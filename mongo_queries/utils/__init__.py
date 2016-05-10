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
  raise RuntimeError("Can't find file creds.json")

def get_creds(user=None) :
  fn = get_creds_fn()
  if os.path.exists(fn) :
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

def get_num_aas_in_chain(pdb_id,chain,db) :
  q = {'pdb_id':pdb_id,'chain_id':chain}
  p = {'resname':1}
  cursor = db.residues_colkeys.find(q,p)
  i = 0
  for res in cursor :
    if res['resname'] in reslist : i += 1
  return i

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
    assert self.has_density_parameters(region=region)
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

  def get_atom_xyz(self, atom) :
    assert isinstance(atom,str)
    if not 'atoms' in self.raw_mongodoc.keys() : return
    atoms = self.raw_mongodoc['atoms']
    if not atom in atoms.keys() : return
    return atoms[atom]["xyz"]

class MongoResidueList(dict) :

  def __init__(self, pdb_id, chain=None) :
    self.pdb_id = pdb_id.lower()
    self.chain = chain

  def get_residues(self, db=None, collection='residues_colkeys') :
    if not db : db = connect(db='pdb_info')
    q = {'pdb_id':self.pdb_id}
    if self.chain : q['chain_id'] = self.chain
    assert hasattr(db,collection)
    dbcol = getattr(db,collection)
    cursor = dbcol.find(q)
    for r in cursor :
      #print r
      self[str(r['_id'])] = MongoResidue(mongodoc = r)
      #if self[r['_id']].resname == 'HOH' : continue
      #if not self[r['_id']].has_density_parameters(): continue
      #print self[r['_id']],self[r['_id']].paases_filter()#,self[r['_id']].next_residue()
      #break

  def get_resolution(self, db=None) :
    if not db : db = connect(db='pdb_info')
    q = {"_id.pdb_id":self.pdb_id.upper()}
    p = {"resolution":1}
    cursor = db.experiment.find(q,p)
    assert cursor.count() == 1
    if "resolution" not in cursor[0].keys() : return
    return cursor[0]["resolution"]

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
