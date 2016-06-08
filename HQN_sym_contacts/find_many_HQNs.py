import os,sys
import argparse
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath,'..'))
from mongo_queries import utils
from nonbonded_sym_contacts import get_nonbonded_sym_contacts
import libtbx
from libtbx import easy_run
import mmtbx.command_line.fetch_pdb
import get_HQN_sym_contacts

_default_distance = 2.5

def parseags() :
  desc = """This script will return a list of many potential inter_asymmetric unit NQH flips from the Top8000.
"""
  parser = argparse.ArgumentParser(description=desc,
                               formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-n', '--number_pdbs', type=int,
                      help="Number of PDBs to analyze")
  s = 'Distance cutoff between the two interacting atoms. default is %.1f'
  parser.add_argument('-d', '--distance_cutoff', type=float,
                       help=s % _default_distance, default=_default_distance)
  args = parser.parse_args()
  return args

def get_pdbs () :
  fn = 'List_hi-res_good-ed.txt'
  pdbs = []
  fle = open(fn,'r')
  for l in fle :
    if l.strip() == '' or l.startswith('#') : continue
    pdbs.append((l[:4],''))
  fle.close()
  return pdbs

def run() :
  args = parseags()
  homology_level = 70
  pdbs = [('1a2z','A')]#utils.get_Top8000_pdb_list(homology_level=homology_level)
  pdbs = get_pdbs()
  s = '%i pdbs found in Top8000 at homology level %i'
  utils.broadcast(s % (len(pdbs),homology_level))
  there = False
  for i,pc in enumerate(pdbs) :
    if i % 100 == 0 : print >> sys.stderr, "Through %i of %i.." % (i,len(pdbs))
    pdb_id,chain = pc
    if pdb_id == '2ykz' :
      there = True
      continue
    if not there : continue
    print pdb_id
    pdbfn = mmtbx.command_line.fetch_pdb.run2(args=[pdb_id])
    assert os.path.exists(pdbfn)
    pairs = get_HQN_sym_contacts.get_HQN_sym_contacts(pdb_file=pdbfn,
                                         distance_cutoff=args.distance_cutoff)
    pairs.write_formatted_pairs()
    # cleanup
    os.remove(pdbfn)
    if args.number_pdbs and i+1 >= args.number_pdbs : break

if __name__ == '__main__' :
  run()
