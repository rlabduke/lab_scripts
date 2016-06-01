import os,sys
import argparse
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath,'../..'))
import mongo_queries.utils

default_resol = 2

def parseags() :
  desc = """This script is meant to identify DNA base pairs between asymmetric units. Give a specific pair to find, e.g. 

$ python get_sym_DNA_bp.py AA 

for adenine : adenine. 

NOTE this script finds nonbonded interactions between given bases which is not gaurenteed to be an authentic base pair."""
  parser = argparse.ArgumentParser(description=desc,
                               formatter_class=argparse.RawTextHelpFormatter)
  s = 'Get pdbs with resolution higher than this. Default is %.1f'
  parser.add_argument('basepair',type=str,help="Specify basepair of interest")
  parser.add_argument('-r','--resolution',type=float,default=default_resol,
                      help= s % default_resol)
  #parser.add_argument('-v','--verbose',action='store_true',help='Be verbose')
  return parser.parse_args()

def run() :
  args = parseags()
  assert len(args.basepair) == 2
  bases = ["A","T","G","C"]
  assert args.basepair[0] in bases and args.basepair[1] in bases,\
            args.basepair[0]
  print args.basepair
 


if __name__ == "__main__" :
  run()
