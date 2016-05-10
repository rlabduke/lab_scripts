import os,sys
import argparse
sys.path.append('..')
import utils

def get_cystines(pdb_id, db) :
  residues = utils.MongoResidueList(pdb_id = pdb_id)
  residues.get_residues(db=db)
  residues.link_residues()
  residue_keys = residues.ordered_keys()
  cystines  = []
  for reskey in residue_keys :
    mongores = residues[reskey]
    if mongores.name == 'CYS' : cystines.append(mongores)
  return cystines

def run(args) :
  desc = "Decribe script here"
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('-r','--resolution',type=float,default=1.2,
                      help='Get pdbs with resolution higher than the given')
  parser.add_argument('-v','--verbose',action='store_true',help='Be verbose')
  args = parser.parse_args()

  # Get pdbs
  pdbs = utils.get_filtered_pdbs(high_resolution=args.resolution)
  s = '%i pdbs found with resolution highr than %.1f'
  utils.broadcast(s % (len(pdbs),args.resolution))

  # Get connection to mongo
  mongocon = utils.MongodbConnection()

  # iterate through pdbs
  for pdb_id in pdbs :
    # Get cystines
    cystines = get_cystines(pdb_id, db=mongocon.db)
    s = '%i cystines found in %s'
    utils.broadcast(s % (len(cystines),pdb_id))
    break

if __name__ == '__main__' :
  run(sys.argv[1:])

