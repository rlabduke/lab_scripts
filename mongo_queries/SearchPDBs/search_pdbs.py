import os,sys
import argparse
sys.path.append('..')
import utils

default_resol = 1.2

def run()
  desc = "This script is meant to identify DNA base pairs between asymmetric "
  desc+= "units. To find a specific pair specify using the -p flag, e.g. "
  desc+= "-p AA for adenine : adenine. Note this script finds nonbonded "
  desc+= "interactions between given bases which is not gaurenteed to be an "
  desc+= "authentic base pair."
  parser = argparse.ArgumentParser(description=desc)
  s = 'Get pdbs with resolution higher than this. Default is %.1f'
  parser.add_argument('-r','--resolution',type=float,default=default_resol,
                      help= s % default_resol)
  parser.add_argument('-v','--verbose',action='store_true',help='Be verbose')
  args = parser.parse_args()

  


if __name__ == "__main__" :
  run()
