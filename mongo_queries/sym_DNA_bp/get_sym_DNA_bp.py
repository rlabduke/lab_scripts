import os,sys
import argparse
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath,'../..'))
import mongo_queries.utils
from nonbonded_sym_contacts import get_nonbonded_sym_contacts

default_resol = 2
default_distance = 2.5

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
  s = 'Distance cutoff between the two interacting atoms. default is %.1f'
  parser.add_argument('-d', '--distance_cutoff', type=float,
                       help=s % default_distance, default=default_distance)
  parser.add_argument('-v','--verbose',action='store_true',help='Be verbose')
  return parser.parse_args()

def get_dna_pdbs(resolution) :
  fn = os.path.join(scriptpath,'dnas.l')
  assert os.path.exists(fn)
  pdbs = []
  fle = open(fn,'r')
  for l in fle :
    if l.strip() == '' or l.startswith('#') : continue
    pdb_id = l.strip()
    assert len(pdb_id) == 4
    if pdb_id not in pdbs : pdbs.append(pdb_id)
  print >> sys.stderr, '\n%i pdbs to iterate.\n' % len(pdbs)
  return pdbs

def fetch_pdb(pdb_id) :
  import mmtbx.command_line.fetch_pdb
  return mmtbx.command_line.fetch_pdb.run2(args=[pdb_id])

def run() :
  args = parseags()
  assert len(args.basepair) == 2
  bases = ["A","T","G","C"]
  assert args.basepair[0] in bases and args.basepair[1] in bases,\
            args.basepair[0]
  filter_residues = ([" D%s" % args.basepair[0]],[" D%s" % args.basepair[1]])
  print >> sys.stderr, "\nFinding %s..." % args.basepair
  # get pdbs with DNA and no RNA
  pdbs = get_dna_pdbs(args.resolution)

  # iterate through pdbs
  thr = "Through %i - %.2f%%"
  found = '%i pair[s] found in %s'
  found_n = 0
  for i,pdb_id in enumerate(pdbs) :
    if i % 100 == 0 : print >> sys.stderr, thr % (i,(i*100.)/len(pdbs))
    pdb_id = '4G0F'
    fn = fetch_pdb(pdb_id)
    if args.verbose : print >> sys.stderr, "Fetched %s" % fn
    reduced_pdb_file = get_nonbonded_sym_contacts.reduce_pdb(fn)
    pairs = get_nonbonded_sym_contacts.NonbondedIinteractions()
    pairs.get_nonbonded_interactions(file_name=reduced_pdb_file,
                                     distance_cutoff=args.distance_cutoff,
                                     filter_residues = filter_residues)
    residue_level_pairs = {}
    for pair in pairs  :
      # We dont care about ribose atoms
      if pair.residues[0].name.strip().endswith('\'') \
         or pair.residues[1].name.strip().endswith('\'') : continue
      key = (pdb_id,pair.residue_ids[0],pair.residue_ids[1])
      if not key in residue_level_pairs.keys() : residue_level_pairs[key] = []
      t = (pair.residues[0].name,pair.residues[1].name,'%.2f'%pair.distance)
      residue_level_pairs[key].append(t) 
    if len(residue_level_pairs) > 0 :
      print >> sys.stderr, found % (len(residue_level_pairs),pdb_id)
      found_n += len(residue_level_pairs)
    for k,l in residue_level_pairs.items() :
      print >> sys.stdout, k
      for e in l : print >> sys.stdout, '    ' + str(e)
    # cleanup
    os.remove(fn)
    os.remove(reduced_pdb_file)
    break
  print >> sys.stderr, '%i TOTAL %s pair[s] found' % (found_n,args.basepair)

if __name__ == "__main__" :
  run()
