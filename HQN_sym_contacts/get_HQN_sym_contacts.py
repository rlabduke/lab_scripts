import os,sys
import argparse
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath,'..'))
from nonbonded_sym_contacts import get_nonbonded_sym_contacts

_default_distance = 2.5
HIS_ATOMS = [" HD2"," HE1"," ND1"," NE2"," HE2"," HD1"]
ASN_ATOMS = ["HD22","HD21"," OD1"]
GLN_ATOMS = ["HE22","HE21"," OE1"]

def parseags() :
  desc = """This script will return all NQH residues contacting anything in a different asymmetric unit.
"""
  parser = argparse.ArgumentParser(description=desc,
                               formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('pdb_file', help='A pdb file')
  s = 'Distance cutoff between the two interacting atoms. default is %.1f'
  parser.add_argument('-d', '--distance_cutoff', type=float,
                       help=s % _default_distance, default=_default_distance)
  args = parser.parse_args()
  assert os.path.exists(args.pdb_file)
  return args

def get_HQN_sym_contacts(pdb_file,distance_cutoff) :
  reduced_pdb_file = get_nonbonded_sym_contacts.reduce_pdb(pdb_file)
  filter_residues = (["HIS","GLN","ASN"],None)
  filter_atoms = ({"HIS":HIS_ATOMS,"ASN":ASN_ATOMS,"GLN":GLN_ATOMS},None)
  pairs = get_nonbonded_sym_contacts.NonbondedIinteractions()
  pairs.get_nonbonded_interactions(file_name=reduced_pdb_file,
                                   distance_cutoff=distance_cutoff,
                                   select = None,
                                   filter_residues=filter_residues,
                                   filter_atoms = filter_atoms)
  return pairs

def run() :
  args = parseags()
  get_HQN_sym_contacts(args.pdb_file, args.distance_cutoff)
  pairs.write_formatted_pairs()

if __name__ == '__main__' :
  run()

