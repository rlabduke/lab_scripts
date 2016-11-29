import os,sys
import argparse
import get_NQH_density
sys.path.append('..')
import utils

def run(args) :
  desc = "Takes a list file of pdb codes (1 code per line) and writes a csv "
  desc+= "for each pdb containing records for each NQH with the following "
  desc+= "headers: pdb_id, chain_id, resseq, icode, altloc, resname, "
  desc+= "worst_rscc, worst_2fo-fc, N, C_O, delta]"
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdb_list', help='A file with pdb codes (one per line)')
  args = parser.parse_args()
  assert os.path.exists(args.pdb_list)

  mongocon = utils.MongodbConnection()
  fle = open(args.pdb_list,'r')
  for l in fle :
    nqhD = get_NQH_density.NQH_density(l.strip())
    nqhD.set_density_values(mongocon.db)
    nfl = open("%s.csv" % l.strip(),'w')
    nqhD.write_csv(out=nfl)
    nfl.close()
    print >> sys.stderr, "%s.csv written." % l.strip()
  fle.close()

if __name__ == '__main__' :
  run(sys.argv[1:])

