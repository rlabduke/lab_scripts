import os,sys
import argparse
sys.path.append('..')
import utils
import numpy

def get_cystines(residues, db) :
  residues.get_residues(db=db)
  residues.link_residues()
  residue_keys = residues.ordered_keys()
  cystines  = []
  for reskey in residue_keys :
    mongores = residues[reskey]
    if mongores.resname == 'CYS' : cystines.append(mongores)
  return cystines

class CysCys(list) :

  def add_cyscys(self,cys0,cys1,resolution,residue_n) :
    self.append({'res0':cys0,'res1':cys1,
                 'resolution':resolution,
                 'residue_n':residue_n})

  def SG_SG_distanced(self,c0,c1) :
    a = numpy.array(c0.get_atom_xyz("SG"))
    b = numpy.array(c1.get_atom_xyz("SG"))
    return numpy.linalg.norm(a-b)

  def write_csv(self,omegatype=None,log=sys.stdout) :
    assert omegatype in [None,'Cis','Trans']
    heads = ['residue0','passes_filter0','residue1','passes_filter1']
    heads+= ['resolution','chain_aa_n','omega','omega_type','SG-SG']
    if len(self) == 0 :
      print >> sys.stderr, "No CysCys here!"
    else :
      print >> log, ','.join(heads)
      for d in self :
        res0 = d['res0'].as_str()
        p0 = d['res0'].passes_filter()
        res1 = d['res1'].as_str()
        p1 = d['res1'].passes_filter()
        #print utils.print_json_pretty(c1.raw_mongodoc)
        if hasattr(d['res1'],'omegalyze') :
          ot = d['res1'].omegalyze.type
          omega = "%.2f" % d['res1'].omegalyze.omega
        else : ot,omega = None,None
        reso = "%.2f" % d['resolution']
        # Get SG-SG distance
        dist = "%.2f" % self.SG_SG_distanced(d['res0'],d['res1'])
        lst = [res0,p0,res1,p1,reso,d['residue_n'],omega,ot,dist]
        assert len(lst) == len(heads)
        print >> log, ','.join([str(e) for e in lst])

def has_ciscys(pdb_id,db) :
  q = {"pdb_id":pdb_id,"resname":"CYS","omegalyze.type":"Cis"}
  cursor = db.residues_colkeys.count(q)
  if cursor == 0 : return False
  return True

def run(args) :
  desc = "Decribe script here"
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('-r','--resolution',type=float,default=1.2,
                      help='Get pdbs with resolution higher than the given')
  parser.add_argument('-v','--verbose',action='store_true',help='Be verbose')
  parser.add_argument('-c','--cisonly',action='store_true',help='Get cis only')
  args = parser.parse_args()

  # Get pdbs
  pdbs = utils.get_filtered_pdbs(high_resolution=args.resolution)
  if args.verbose :
    s = '%i pdbs found with resolution highr than %.1f'
    utils.broadcast(s % (len(pdbs),args.resolution))

  # Get connection to mongo
  mongocon = utils.MongodbConnection()

  pdbs = ['3edh','2z63','2z66','3zui','4mge','3sr3','2q3z','3dst',
          '3dsu','1wd3','3h0t','3h0u','4nn5','4nf4','3hol','3t6q']
  # iterate through pdbs
  cc = CysCys()
  for i,pdb_id in enumerate(pdbs) :
    if i % 100 == 0 : print >> sys.stderr, "Through %i of %i.." % (i,len(pdbs))
    pdb_id = pdb_id.lower()#'3hol'
    if args.cisonly :
      if not has_ciscys(pdb_id,db=mongocon.db) : continue
    #print has_ciscys(pdb_id,db=mongocon.db)
    # Get cystines
    residues = utils.MongoResidueList(pdb_id = pdb_id)
    cystines = get_cystines(residues, db=mongocon.db)
    if args.verbose :
      s = '%i cystines found in %s'
      utils.broadcast(s % (len(cystines),pdb_id))
    # See if any cystine is preceeded by cystine
    for cys in cystines :
      # Get previous residues if they are cystine
      precedes = [res for res in cys.prevres if res.resname == 'CYS' ]
      if len(precedes) == 0 : continue
      # get resolution
      resolution = residues.get_resolution(db=mongocon.db)
      # get number of residues
      residue_n = utils.get_num_aas_in_chain(pdb_id,cys.chain_id,db=mongocon.db)
      # add the cyscys
      for res in precedes : cc.add_cyscys(res,cys,resolution,residue_n)
    if args.verbose :
      s = '%i adjacent cystines found in %s'
      utils.broadcast(s % (len(precedes),pdb_id))
    #if len(cc) == 30 : break
    #break
  cc.write_csv()

if __name__ == '__main__' :
  run(sys.argv[1:])

