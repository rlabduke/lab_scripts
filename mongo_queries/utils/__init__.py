import os,sys
from pymongo import MongoClient
import json
import re
import pprint
import math
import time

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

  def __init__(self,db_name='pdb_info',user='rlaber') :
    self.user, self.pwd = get_creds(user)
    print >> sys.stderr, 'Logged into Mongo as %s.' % self.user
    self.db_name = db_name
    self.connect()

  def connect(self) :
    uri = "mongodb://%s:%s@127.0.0.1/?authSource=test"
    self.client = MongoClient(uri % (self.user,self.pwd))
    print 
    assert hasattr(self.client,self.db_name),'%s not in MongoDB'%self.db_name
    self.set_db(self.db_name)

  def set_db(self,db) :
    self.db_name = db
    self.db = getattr(self.client,db)

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
                      verbose=False,
                      low_resolution=0.0) :
  assert type(high_resolution) == float
  if connection : isinstance(connection,MongodbConnection)
  else : connection = MongodbConnection()
  qd = {"experimental_method":experimental_method}
  qd['resolution'] = {'$lte':high_resolution, '$gte':low_resolution}
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

def get_Top8000_pdb_list(homology_level=70,verbose=False,connection=None) :
  # Although not required, it is generally a good idea to provida a connection,
  # which is a MongodbConnection object, as you will likely use this object
  # later in your script as you iterate through the pdbs gotten here.
  assert homology_level in [50,70,90,95]
  q = {"in_mtz_%i"%homology_level:1}
  # we only want pdb_id and chain
  p = {"pdb_id":1,"chain":1}
  if not connection : mongocon = MongodbConnection()
  else : mongocon = connection
  mongocon.set_db(db='top8000_rota_data')
  cursor = mongocon.db.versions_2.find(q,p)
  if verbose : print >> sys.stderr, 'fetched %i PDB-chains.' % cursor.count()
  return [(e["pdb_id"],e["chain"]) for e in cursor]

hetero_list = []
hydrogen_list = []
#print(os.path.join(os.path.dirname(os.path.realpath(__file__)),"sc-connect.props"))
with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),"sc-connect.props")) as pdb_atoms_file:
  hetero_list = []
  hydrogen_list = []
  for line in pdb_atoms_file:
    if not line.startswith("#") and "=" in line:
      split_line = line.split(" = ")
      split_line[1]=split_line[1].strip()
      split_line[1]=split_line[1].strip("\"")
      split_line[1]=split_line[1].strip(";")
      #print split_line
      if ".hy" in split_line[0]:
        atom_pairs = split_line[1].split(";")
        for atom_pair in atom_pairs:
          atoms = atom_pair.split(",")
          if not atoms[1] in hydrogen_list:
            hydrogen_list.append(atoms[1])
      else:
        atoms = re.split(',|;', split_line[1])
        for atom in atoms:
          if not atom in hetero_list:
            hetero_list.append(atom)
  #print hydrogen_list
  hetero_list = hetero_list + hydrogen_list
  #print hetero_list
  
def atom_sort(atom):
  atom = format_atom(atom)
  if atom in hetero_list:
    return hetero_list.index(atom)
  return -1
  
def format_atom(atom):
  if len(atom) == 4: return atom
  if len(atom) == 1 or len(atom) == 2: return "{:^4}".format(atom)
  if len(atom) == 3: return "{:>4}".format(atom)
  return "????"

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
    self.translated_atoms = self.init_translated_coords() # should be a deep copy
    self.set_filter_thresholds()
    self.set_omega()
    self.set_rama()

  def __str__(self) :
    return self.as_str()
    
  def clone(self):
    cls = self.__class__
    residue_clone = cls.__new__(cls)
    residue_clone.__dict__.update(self.__dict__)
    residue_clone.translated_atoms = self.init_translated_coords()
    return residue_clone

  def as_str(self,noalt=False) :
    if noalt :
      return ' '.join([self.pdb_id,self.model_id,self.chain_id,
                       self.resseq,self.icode,self.resname])
    return ' '.join([self.pdb_id,self.model_id,self.chain_id,
                     self.altloc,self.resseq,self.icode,self.resname])

  def init_translated_coords(self):
    if not 'atoms' in self.raw_mongodoc.keys() : return
    trans_atoms = {}
    atoms = self.raw_mongodoc['atoms']
    for atom in atoms:
      atom_dict_copy = {}
      xyz_copy = []
      for key in atoms[atom]:
        if key == 'xyz':
          xyz = atoms[atom]['xyz']
          for coord in xyz:
            xyz_copy.append(coord)
            atom_dict_copy['xyz'] = xyz_copy
        else:
          atom_dict_copy[key] = atoms[atom][key]
      trans_atoms[atom] = atom_dict_copy
    return trans_atoms

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
    
  def is_outlier(self):
    analysis_types = ['rotalyze', 'omegalyze', 'ramalyze']
    mongo_keys = self.raw_mongodoc.keys()
    #pprint.pprint(mongo_keys)
    if not 'ramalyze' in mongo_keys: 
      return True # knock out chain ends
    if 'clashes' in mongo_keys: 
      #print "clashes found"
      return True
    if "angle" in mongo_keys: 
      #print "angleoutfound"
      return True
    if "bondlength" in mongo_keys: 
      #print "bondoutfound"
      return True
    for analysis in analysis_types:
      if analysis in mongo_keys:
        if self.raw_mongodoc[analysis]['is_outlier']: return True
    if "worst_bb" in mongo_keys:
      #print "worst_bbfound"
      try:
        if self.raw_mongodoc['worst_bb']['adp']['value'] > 30: return True
      except KeyError:
        print(self.pdb_id+self.resname+self.resseq+" seems to be missing an adp value")
    if "cablam" in mongo_keys:
      #print "cablam found"
      if self.raw_mongodoc['cablam'][0]['c_alpha_geom_outlier']: return True
    bb_atoms = ["N", "C", "CA", "O"]
    atoms = self.raw_mongodoc['atoms']
    for bb_atom in bb_atoms:
      if not bb_atom in atoms:
        return True
    

  def set_omega(self) :
    if 'omegalyze' not in self.raw_mongodoc.keys() : return
    self.omegalyze = group_args(**self.raw_mongodoc['omegalyze'])

  def set_rama(self) :
    if 'ramalyze' not in self.raw_mongodoc.keys() : return
    self.ramalyze = group_args(**self.raw_mongodoc['ramalyze'])

  def get_atoms(self):
    if not 'atoms' in self.raw_mongodoc.keys() : return
    atoms = self.raw_mongodoc['atoms']
    return atoms

  def get_atom_xyz(self, atom) :
    assert isinstance(atom,str)
    if not 'atoms' in self.raw_mongodoc.keys() : return
    atoms = self.raw_mongodoc['atoms']
    if not atom in atoms.keys() : return
    return atoms[atom]["xyz"]
    
  def get_translated_xyz(self, atom):
    if not atom in self.translated_atoms: return
    return self.translated_atoms[atom]["xyz"]
    
  def set_atom_xyz(self, atom, xyz):
    assert isinstance(atom, str)
    #if not 'atoms' in self.raw_mongodoc.keys() : return
    #print("pre set translated_atoms")
    #pprint.pprint(self.translated_atoms)
    if not atom in self.translated_atoms.keys() : return
    #atoms_str = str(self.raw_mongodoc['atoms'])
    #print("setting translated atoms")
    self.translated_atoms[atom]["xyz"] = list(xyz)
    #print("post set translated_atoms")
    #pprint.pprint(self.translated_atoms)
    #assert atoms_str == str(self.raw_mongodoc['atoms']) # test to make sure original atoms not changed
        
  def lowest_occ(self) :
    lowest_occ = 1
    atoms = self.raw_mongodoc['atoms']
    for an,d in atoms.items() :
      if d['occ'] < lowest_occ : lowest_occ = d['occ']
    return lowest_occ
  
  # attempt to figure out the element type from the atom name
  def get_atom_element(self, atom):
    if len(atom) ==1: return atom
    if "C" in atom: return "C"
    if "H" in atom: return "H"
    if "N" in atom: return "N"
    if "O" in atom: return "O"
    if "S" in atom: return "S"
    
  # reconstruct the residue's atom records?  
  def get_atom_record(self, atom, atoms, atom_number, renumber_resnum=None):
    assert isinstance(atom,str)
    #if not 'atoms' in self.raw_mongodoc.keys() : return
    #atoms = self.raw_mongodoc['atoms']
    if not atom in atoms.keys() : return
    if renumber_resnum is None:
      renumber_resnum = self.resseq
    atom_formatter = "ATOM  {atomnum:>5} {atomname}{altloc:>1}{resname} {chainid:>1}{resnum:>4}{icode:>1}   {xcoord:>8.3f}{ycoord:>8.3f}{zcoord:>8.3f}{occ:>6.2f}{bfact:>6.2f}          {element:>2}  {extra}"
    return atom_formatter.format(atomnum=atom_number,
                                 atomname=format_atom(atom),
                                 altloc=self.altloc,
                                 resname=self.resname,
                                 chainid=self.chain_id,
                                 icode=self.icode,
                                 resnum=renumber_resnum,
                                 xcoord=atoms[atom]["xyz"][0],
                                 ycoord=atoms[atom]["xyz"][1],
                                 zcoord=atoms[atom]["xyz"][2],
                                 occ=atoms[atom]["occ"],
                                 bfact=atoms[atom]["adp"],
                                 element=self.get_atom_element(atom),
                                 extra=self.pdb_id+self.resname+self.resseq)
                                     
  #translated_atoms must have actual atom names
  def get_atom_records(self, translated=False, region="all", renumber_residue=None):
    assert region in ['all','bb']
    bb_atom = ["N", "C", "CA", "O", "CB"]
    if not 'atoms' in self.raw_mongodoc.keys() : return
    atoms = {}
    if translated:
      atoms = self.translated_atoms
    else:
      atoms = self.raw_mongodoc['atoms']
    all_records=""
    for atom in sorted(atoms.keys(), key=atom_sort):
      if region=='bb' and atom in bb_atom or region=='all':
        all_records = all_records+(self.get_atom_record(str(atom), atoms, "1", renumber_residue)+"\n")
    return all_records
    
# this is mainly intended to be a list of the residues in a fragment of a PDB
# for use for generating tom's fragments
# HACK: only works for fragments up to 9 in length
class MongoPdbFragment(object):
  
  def __init__(self, mongo_residues):
    self.residues = mongo_residues
    bb_atom = ["N", "C", "CA", "O", "CB"]
    bb_atoms_dict = {}
    for i, mongo_res in enumerate(mongo_residues):
      for atom in bb_atom:
        xyz = mongo_res.get_atom_xyz(atom)
        if atom != "CB" or not xyz is None:
          bb_atoms_dict[str(i)+atom] = list(xyz)
    self.bb_atom_coords = bb_atoms_dict
    
  def __str__(self):
    return str(self.residues[0])+" - "+str(self.residues[-1])
  
  def __repr__(self):
    return str(self.residues[0])+" - "+str(self.residues[-1])
    
  # bb_atom_coords has the format: {"0N":[x,y,z], "0CA":[x,y,z], ......} 
  def get_bb_atoms(self):
    return self.bb_atom_coords
    
  def set_bb_atoms(self, bb_atoms_dict):
    self.bb_atom_coords=bb_atoms_dict
    for i, mongo_res in enumerate(self.residues):
      #print(mongo_res)
      for bb_atom in ["N", "C", "CA", "O", "CB"]:
        #print(atom[1:])
        if str(i)+bb_atom in bb_atoms_dict:
          mongo_res.set_atom_xyz(bb_atom,bb_atoms_dict[str(i)+bb_atom])
        #print("translatedatom: "+str(mongo_res.get_translated_xyz(bb_atom)))
        
  def get_atom_records(self, translated=False, region="all"):
    records=""
    for i, residue in enumerate(self.residues):
      records = records+residue.get_atom_records(translated, region, i+1)
    return records
    
  def get_residues(self):
    return self.residues
    
  def get_rmsd(self, fragment):
    assert len(self.residues) == len(fragment.get_residues())
    ref_fragment_coords = self.bb_atom_coords
    test_fragment_coords = fragment.get_bb_atoms()
    atom_counter = 0
    dist_sum = 0
    for i, mongo_res in enumerate(self.residues):
      for bb_atom in ["N", "C", "CA", "O"]:
        current_atom_key = str(i)+bb_atom
        if current_atom_key in ref_fragment_coords:
          if current_atom_key in test_fragment_coords:
            atom_counter = atom_counter + 1
            ref_xyz = ref_fragment_coords[current_atom_key]
            test_xyz = test_fragment_coords[current_atom_key]
            #dist_sum = dist_sum + ((ref_xyz[0]-test_xyz[0])**2+(ref_xyz[1]-test_xyz[1])**2+(ref_xyz[2]-test_xyz[2])**2)
            dist_sum = dist_sum + (math.pow((ref_xyz[0]-test_xyz[0]),2)+math.pow((ref_xyz[1]-test_xyz[1]),2)+math.pow((ref_xyz[2]-test_xyz[2]),2))
    rmsd = math.sqrt(dist_sum/atom_counter)
    #print("rmsd: "+str(rmsd))
    return rmsd

class MongoResidueList(dict) :

  def __init__(self, db, pdb_id, chain=None, collection='residues_colkeys') :
    self.pdb_id = pdb_id.lower()
    self.chain = chain
    self.db = db
    #start_time = time.time()
    self.get_residues(collection)
    #elapsed_time = time.time() - start_time
    #print str(elapsed_time) + " time taken to get residues"
    self.link_residues()
    #self.set_counts()

  def get_residues(self, collection='residues_colkeys') :
    q = {'pdb_id':self.pdb_id}
    if self.chain : q['chain_id'] = self.chain
    #q['clashes']={"$exists":False}
    #q['angle']={"$exists":False}
    #q['bondlength']={"$exists":False}
    #q['ramalyze']={"$exists":True}
    assert hasattr(self.db,collection)
    dbcol = getattr(self.db,collection)
    cursor = dbcol.find(q)
    #cursor = dbcol.find(q).explain()
    #pprint.pprint(cursor)
    #print cursor.count()
    #cursor2 = dbcol.find_one()
    #print cursor2
    for r in cursor :
      #print r
      #start_time = time.time()
      self[str(r['_id'])] = MongoResidue(mongodoc = r)
      #elapsed_time = time.time() - start_time
      #print str(elapsed_time) + " time taken to loop residues"
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
