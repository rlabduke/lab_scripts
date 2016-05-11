import os,sys
import argparse
sys.path.append('..')
import utils
import numpy

class NonProCis(list) :

  def add_nonprocis(self) : pass

def has_non_pro_cis(pdb_id,chain,db) :
  q = {"pdb_id":pdb_id,"chain":chain,"resname":{"$ne":"PRO"}}
  q["omegalyze.type"] = "CIS"
  if db.residues_colkeys.count(q) == 0 : return False
  return True

def get_npc_residues(residues) :
  residue_keys = residues.ordered_keys()
  npc_residues  = []
  for reskey in residue_keys :
    mongores = residues[reskey]
    if not hasattr(mongores,'omegalyze') : continue
    if mongores.omegalyze.type == "Cis" and mongores.resname != 'CYS' :
      npc_residues.append(mongores)
  return npc_residues

def run(args) :
  desc = "A query script to getnon-pro cis residues from the Top8000 at a given"
  desc+= " homology level."
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('-o','--homology_level', type=int,default=70,
                   help='Homology level can be 50, 70, 90, or 95. Default=70')
  parser.add_argument('-v','--verbose',action='store_true',help='Be verbose')
  args = parser.parse_args()

  # Get pdbs
  pdbs = utils.get_Top8000_pdb_list(homology_level=args.homology_level)
  if args.verbose :
    s = '%i pdbs found in Top8000 at homology level %i'
    utils.broadcast(s % (len(pdbs),args.homology_level))

  # Get connection to mongo
  mongocon = utils.MongodbConnection()

  # Iterate through pdbs
  counts = {"all_aa":{'n_unique':0,'n_alts':0,'passfilter':0,'alt_passfilter':0},
          "pro":{'n_unique':0,'n_alts':0,'passfilter':0,'alt_passfilter':0},
          "cispro":{'n_unique':0,'n_alts':0,'passfilter':0,'alt_passfilter':0},
          "nonpro":{'n_unique':0,'n_alts':0,'passfilter':0,'alt_passfilter':0},
          "nonprocis":{'n_unique':0,'n_alts':0,'passfilter':0,'alt_passfilter':0}}
  nonprocis = NonProCis()
  for i,pc in enumerate(pdbs) :
    if i % 100 == 0 : print >> sys.stderr, "Through %i of %i.." % (i,len(pdbs))
    pdb_id,chain = pc
    pdb_id,chain = "193l","A"
    # get number of canonical aas in given pdb_chain
    counts['all_aa']['n_unique'] += utils.get_num_aas_in_chain(pdb_id,chain,mongocon.db)
    print counts['all_aa']['n_unique']
    residues = utils.MongoResidueList(db=mongocon.db, pdb_id=pdb_id, chain=chain)
    print residues.counts.unique_canonical_aa
    print residues.counts.unique_canonical_aa_filter
    print residues.counts.all_canonical_aa
    print residues.counts.all_canonical_aa_filter
    break
    # skip if there are no non-pro cis in pdb_id,chain
    if not has_non_pro_cis(pdb_id,chain,mongocon.db) : continue
    # get non-pro cis residues
    residues = utils.MongoResidueList(db=mongocon.db,pdb_id=pdb_id, chain=chain)
    npc_residues = get_npc_residues(residues)
    # Sanity check
    assert len(npc_residues) > 0
    
    if i > 5 : break

if __name__ == '__main__' :
  run(sys.argv[1:])

