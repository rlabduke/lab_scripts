import os,sys
from pymongo import MongoClient
import json

reslist = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 
           'PHE', 'ILE', 'HIS', 'LYS', 'MET',
           'LEU', 'ASN', 'GLN', 'PRO', 'SER',
           'ARG', 'THR', 'TRP', 'VAL', 'TYR', 'MSE']

collections = ['pdb_residues','file_info','experiment']

def print_json_pretty(d,log=sys.stdout) :
  try :
    print >> log, json.dumps(d,indent=4, separators=(',', ': '))
  except :
    print >> log, 'Could not print json document'

def broadcast(message,log=sys.stderr) :
  print >> log, "*"*79 + '\n'
  print >> log, "%s" % message
  print >> log, "\n" + "*"*79

def broadcast_query(query,collection,n=None,project=None,log=sys.stderr) :
  assert type(query) == dict
  print >> log, '*' * 30 + '  query summary  ' + '*' * 30
  print >> log, 'Collection : %s' % collection
  print >> log, 'Query :'
  print_json_pretty(query)
  if project :
    assert type(project) == dict
    print >> log, 'Project :'
    print_json_pretty(project)
  if n :
    assert type(n) == int
    print >> log, 'Number of documents returned : %i' % n

def get_creds_fn() :
  for dirname in sys.path:
    candidate = os.path.join(dirname, 'creds.json')
    if os.path.isfile(candidate):
      return candidate
  #raise RuntimeError("Can't find file creds.json")
  return None

def get_creds(user=None) :
  fn = get_creds_fn()
  if fn and os.path.exists(fn) :
    fle=open(fn,'r')
    d=json.load(fle)
    fle.close()
    return d['user'],d['dwp']
  import getpass
  if not user :
    user = getpass.getuser()
  print >> sys.stderr, "Please enter password for %s :" % user
  pwd = getpass.getpass()
  return user,pwd

class group_args(object):

  def __init__(self, **keyword_arguments):
    self.__dict__.update(keyword_arguments)

  def __call__(self):
    return self.__dict__

class MongodbConnection(object) :

  def __init__(self,db_name='pdb_info',user=None) :
    self.user, self.pwd = get_creds(user)
    self.db_name = db_name
    self.connect()

  def connect(self) :
    uri = "mongodb://%s:%s@daneel.research.duhs.duke.edu/"
    client = MongoClient(uri % (self.user,self.pwd))
    assert hasattr(client,self.db_name), '%s not on daneel'  % self.db_name
    self.db = getattr(client,self.db_name)

def get_num_aas_in_chain(pdb_id,chain,db,return_both_unique_alt=False) :
  q = {'pdb_id':pdb_id,'chain_id':chain}
  p = {"chain_id":1,"resseq":1,'resname':1,"icode":1,'altloc':1}
  cursor = db.residues_colkeys.find(q,p)
  i = 0
  resl = []
  for d in cursor :
    if d['resname'] in reslist :
      i += 1
      s = ''.join([d["chain_id"],str(d["resseq"]),d['resname'],d["icode"]])
      if s not in resl : resl.append(s)
  if return_both_unique_alt : return len(resl),i
  return len(resl)

def get_filtered_pdbs(high_resolution,
                      connection=None,
                      limit=None,
                      experimental_method="X-ray diffraction",
                      verbose=False) :
  assert type(high_resolution) == float
  if connection : isinstance(connection,MongodbConnection)
  else : connection = MongodbConnection()
  qd = {"experimental_method":experimental_method}
  qd['resolution'] = {'$lte':high_resolution}
  pd = {"_id":1}
  if limit :
    assert type(limit) == int
    cursor = connection.db.experiment.find(qd,pd).limit(limit)
  else : cursor = connection.db.experiment.find(qd,pd)
  pdbs = []
  for d in cursor :
    pdbs.append(d["_id"]["pdb_id"])
  if verbose :
    broadcast_query(query= qd,collection='experiment',n=len(pdbs),project=pd)
  return pdbs

def get_Top8000_pdb_list(homology_level=70,verbose=False) :
  assert homology_level in [50,70,90,95]
  q = {"in_mtz_%i"%homology_level:1}
  # we only want pdb_id and chain
  p = {"pdb_id":1,"chain":1}
  mongocon = MongodbConnection(db_name='top8000_rota_data')
  cursor = mongocon.db.versions_2.find(q,p)
  if verbose : print >> sys.stderr, 'fetched %i PDB-chains.' % cursor.count()
  return [(e["pdb_id"],e["chain"]) for e in cursor]

class MongoResidue(object) :

  def __init__(self, mongodoc) :
    self.raw_mongodoc = mongodoc
    self.id = self.raw_mongodoc['_id']
    sstr = self.id.split(':')
    assert len(sstr) == 7
    self.pdb_id   = sstr[0]
    self.model_id = sstr[1]
    self.chain_id = sstr[2]
    self.icode    = sstr[3]
    self.resseq   = sstr[4]
    self.altloc   = sstr[5]
    self.resname  = sstr[6]
    self.nextres  = []
    self.prevres  = []
    self.set_filter_thresholds()
    self.set_omega()
    self.set_rama()

  def __str__(self) :
    return self.as_str()

  def as_str(self,noalt=False) :
    if noalt :
      return ' '.join([self.pdb_id,self.model_id,self.chain_id,
                       self.resseq,self.icode,self.resname])
    return ' '.join([self.pdb_id,self.model_id,self.chain_id,
                     self.altloc,self.resseq,self.icode,self.resname])

  def next_residue_key_list(self) :
    if 'cablam' not in self.raw_mongodoc.keys() : return
    nrs = [cablam['next'] for cablam in self.raw_mongodoc['cablam']]
    return nrs

  def prev_residue_key_list(self) :
    if 'cablam' not in self.raw_mongodoc.keys() : return
    nrs = [cablam['prev'] for cablam in self.raw_mongodoc['cablam']]
    return nrs

  def deposit_next_residue(self,nextres) :
    if not nextres in self.nextres : self.nextres.append(nextres)

  def deposit_prev_residue(self,prevres) :
    if not prevres in self.prevres : self.prevres.append(prevres)

  def set_filter_thresholds(self,rscc=0.7,mapvalue=1.1,adp=40) :
    self.rscc_threshold     = rscc
    self.mapvalue_threshold = mapvalue
    self.adp_threshold      = adp

  def has_density_parameters(self,region='all') :
    return 'worst_%s' % region in self.raw_mongodoc.keys()

  def passes_filter(self,region='all') :
    assert region in ['all','sc','bb']
    worstregion = 'worst_%s' % region
    #assert self.has_density_parameters(region=region), self.as_str()
    if not self.has_density_parameters(region=region) : return False
    if 'twoFo_DFc_value' not in self.raw_mongodoc[worstregion] : return False
    if self.raw_mongodoc[worstregion]['twoFo_DFc_value']['value'] < \
            self.mapvalue_threshold : return False
    if self.raw_mongodoc[worstregion]['adp']['value'] > \
            self.adp_threshold : return False
    if self.raw_mongodoc[worstregion]['rscc']['value'] < \
            self.rscc_threshold : return False
    return True

  def set_omega(self) :
    if 'omegalyze' not in self.raw_mongodoc.keys() : return
    self.omegalyze = group_args(**self.raw_mongodoc['omegalyze'])

  def set_rama(self) :
    if 'ramalyze' not in self.raw_mongodoc.keys() : return
    self.ramalyze = group_args(**self.raw_mongodoc['ramalyze'])

  def get_atom_xyz(self, atom) :
    assert isinstance(atom,str)
    if not 'atoms' in self.raw_mongodoc.keys() : return
    atoms = self.raw_mongodoc['atoms']
    if not atom in atoms.keys() : return
    return atoms[atom]["xyz"]

  def lowest_occ(self) :
    lowest_occ = 1
    atoms = self.raw_mongodoc['atoms']
    for an,d in atoms.items() :
      if d['occ'] < lowest_occ : lowest_occ = d['occ']
    return lowest_occ

class MongoResidueList(dict) :

  def __init__(self, db, pdb_id, chain=None) :
    self.pdb_id = pdb_id.lower()
    self.chain = chain
    self.db = db
    self.get_residues()
    self.link_residues()
    #self.set_counts()

  def get_residues(self, collection='residues_colkeys') :
    q = {'pdb_id':self.pdb_id}
    if self.chain : q['chain_id'] = self.chain
    assert hasattr(self.db,collection)
    dbcol = getattr(self.db,collection)
    cursor = dbcol.find(q)
    for r in cursor :
      #print r
      self[str(r['_id'])] = MongoResidue(mongodoc = r)
      #if self[r['_id']].resname == 'HOH' : continue
      #if not self[r['_id']].has_density_parameters(): continue
      #print self[r['_id']],self[r['_id']].paases_filter()#,self[r['_id']].next_residue()
      #break

  def get_resolution(self) :
    q = {"_id.pdb_id":self.pdb_id.upper()}
    p = {"resolution":1}
    cursor = self.db.experiment.find(q,p)
    if cursor.count() == 0 : return 0
    if "resolution" not in cursor[0].keys() : return
    return cursor[0]["resolution"]

  def set_counts(self) :
    self.counts = group_args(unique_canonical_aa = 0,
                            all_canonical_aa = 0,
                            unique_canonical_aa_filter = 0,
                            all_canonical_aa_filter = 0)
    keys = self.ordered_keys()
    for k in keys :
       mongores = self[k]
       noalt = mongores.as_str(noalt=True)
       if mongores.resname not in reslist : continue
       if mongores.passes_filter(region='bb') : passes = True
       else : passes = False
       self.counts.all_canonical_aa += 1
       if passes : self.counts.all_canonical_aa_filter += 1
       if noalt not in canon :
         self.counts.unique_canonical_aa += 1
       if passes : 
         self.counts.unique_canonical_aa_filter += 1*mongores.lowest_occ()

  def ordered_keys(self) :
    keys = self.keys()
    keys.sort()
    return keys

  def drop_alt(self,reskey) :
    sk = reskey.split(':')
    sk[5] = ''
    return ':'.join(sk)

  def get_alt_keys(self,reskey) :
    keys = []
    sk = reskey.split(':')
    for s in '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz' :
      nl = sk[:]
      nl[5] = s
      keys.append(':'.join(nl))
    return keys

  def link_residues(self) :
    keys = self.ordered_keys()
    for k in keys :
      prev_residues = self[k].prev_residue_key_list()
      next_residues = self[k].next_residue_key_list()
      if not prev_residues and not next_residues : continue
      #print k, self[k].prev_residue_key_list(), self[k].next_residue_key_list()
      if prev_residues :
        for reskey in prev_residues :
          if not reskey : continue
          if reskey in self.keys() : self[k].deposit_prev_residue(self[reskey])
          elif reskey.split(':')[5] != '' :
            # alt exists, probably a cablam alt and not a model alt
            newkey = self.drop_alt(reskey)
            if newkey in self.keys(): self[k].deposit_prev_residue(self[newkey])
            else :
              print >> sys.stderr, 'Cannot find previous key noalt %s' % reskey
          else :
            # the next residue has alts. This occurs when a residue only has
            # alts in sc thus there are no cablam (where we get adjacency) alts.
            newkeys = self.get_alt_keys(reskey)
            added = 0
            for nk in newkeys :
              if nk in self.keys():
                self[k].deposit_prev_residue(self[nk])
                added += 1
            if added == 0 :
              print >> sys.stderr, 'Cannot find previous key alt %s' % reskey
      if next_residues :
        for reskey in next_residues :
          if not reskey : continue
          if reskey in self.keys() : self[k].deposit_next_residue(self[reskey])
          elif reskey.split(':')[5] != '' :
            # alt exists, probably a cablam alt amd not a model alt
            newkey = self.drop_alt(reskey)
            if newkey in self.keys(): self[k].deposit_next_residue(self[newkey])
            else :
              print >> sys.stderr, 'Cannot find next key noalt %s' % reskey
          else :
            # the next residue has alts. This occurs when a residue only has
            # alts in sc thus there are no cablam (where we get adjacency) alts.
            newkeys = self.get_alt_keys(reskey)
            added = 0
            for nk in newkeys :
              if nk in self.keys():
                self[k].deposit_next_residue(self[nk])
                added += 1
            if added == 0 :
              print >> sys.stderr, 'Cannot find next key alt %s' % reskey
