import os,sys
import argparse
sys.path.append('..')
from utils import pdbe_utils
import mongo_utils

def get_residues_colkeys_pdbs(mcon,limit=None) :
  pdbs = []
  cursor = mcon.db.residues_colkeys.distinct("pdb_id",{})
  for i,_id in enumerate(cursor) :
    #print _id
    pdbs.append(_id)
  return pdbs

def run(args) :
  desc = "This script takes a collection in the pdb_info database on caldera "
  desc+= "and print the missing pdb ids (one per line) to stdout. The vision "
  desc+= "is for this script to be imported and there functiond herein to be "
  desc+= "used in other scripts."
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('collection', help='A collection in pdb_info on caldera')
  parser.add_argument('-o', dest='outfile' ,help='Output file name.')
  args = parser.parse_args()

  current_collections = ["residues_colkeys"]
  assert args.collection in current_collections

  # connect to mongo on caldera
  mcon = mongo_utils.MongodbConnection()

  if args.collection == "residues_colkeys" :
    pdbs = get_residues_colkeys_pdbs(mcon)
    if args.outfile : fle = open(args.outfile,'w')
    else : fle = sys.stdout
    for i,e in enumerate(pdbs) :
      print >> fle, e
      if not args.outfile and i == 10 :
        # If not writing to an outfile, don't write all pdb codes
        print "...%i more pdbs" % (len(pdbs)-i)
        break
    if args.outfile :
      fle.close()
      print >> sys.stderr, '\n%s written.\n' % args.outfile

if __name__ == '__main__' :
  run(sys.argv[1:])

