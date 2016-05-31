import os,sys
import argparse
import libtbx.load_env
from mmtbx import monomer_library
from mmtbx.monomer_library import pdb_interpretation
from libtbx import group_args
from libtbx import easy_run
import re

def parse_label(lab) :
  assert isinstance(lab,str),type(lab)
  ab = re.compile('pdb=".{15}"')
  assert ab.match(lab),lab
  lab = lab[5:-1]
  return group_args(
    name = lab[:4],
    alt_loc = lab[4],
    resname = lab[5:8],
    chain = lab[8:10],
    resseq = lab[10:14],
    icode = lab[14],
    id = lab)

class nonbondedPair(object) :

  def __init__(self,nonb,scatterers) :
    labels,i_seq,j_seq,self.distance,\
       self.vdw_distance,self.sym_op_j,self.rt_mx=nonb
    self.scatterers = scatterers
    self.set_pair_labels(i_seq,j_seq)

  def __str__(self) :
    return "%s - %.2f - %s"% \
        (self.residues[0].id,self.distance,self.residues[1].id)

  def set_pair_labels(self,i_seq,j_seq) :
    ilab = parse_label(self.scatterers[i_seq].label)
    jlab = parse_label(self.scatterers[j_seq].label)
    self.residues = (ilab,jlab)

def reduce_pdb(pdb_file) :
  # returns the new filename
  assert pdb_file.endswith(".pdb") #TODO: Can we do CIF files?
  assert libtbx.env.has_module(name="reduce")
  cmd = "phenix.reduce %s" % pdb_file
  out = easy_run.fully_buffered(cmd)
  assert out.return_code == 0
  nfn = pdb_file.replace(".pdb","H.pdb")
  fle = open(nfn,'w')
  for l in out.stdout_lines: print >> fle, l
  fle.close()
  print >> sys.stderr, '%s written...' % nfn
  return nfn

def get_nonbonded_interactions(file_name, distance_cutoff = 2) :
  # get geostd
  chem_data = libtbx.env.find_in_repositories(
    relative_path="chem_data/geostd",
    test=os.path.isdir)
  # setup monomer_library server
  if chem_data is not None:
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
  # processed pdb file
  pdb_processed_file = pdb_interpretation.process(
    file_name=file_name,
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    force_symmetry = True)
  # get grm
  grm = pdb_processed_file.geometry_restraints_manager(
    show_energies = False,
    plain_pairs_radius = 5.0)
  xrs = pdb_processed_file.xray_structure()
  sites_cart = xrs.sites_cart()
  nonbonded_proxies = grm.pair_proxies().nonbonded_proxies
  get_sorted_result = nonbonded_proxies.get_sorted(
        by_value="delta",
        sites_cart=sites_cart)
  if get_sorted_result is None : return
  sorted_nonb, n_not_shown = get_sorted_result
  # iterate through nonbonded pairs and skip ones < distance and ones without 
  # symmetry
  n_nonb = len(sorted_nonb)
  i = 0
  nonb_pairs = []
  while i < n_nonb and sorted_nonb[i][3] < distance_cutoff :
    if sorted_nonb[i][-1] is None : i += 1; continue
    nonb_pairs.append(nonbondedPair(sorted_nonb[i],xrs.scatterers()))
    i += 1
  return nonb_pairs

def run(args) :
  desc = "This script will show you inter-asymmetric-unit contacts from a PDB "
  desc+= "with a CRYST card. I suggest 3g6t as it has a nice base pair between "
  desc+= " asymmetric units, C 1 to C 2. This script adds hysrogens for you. "
  desc+= "\n\nAGAIN, THESE ARE ONLY INTER_ASYMMETRIC_UNIT INTERACTIONS!!\n\n"
  parser = argparse.ArgumentParser(description=desc,
                                  formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('pdb_file', help='A pdb file')
  s = 'Distance cutoff between the two interacting atoms'
  parser.add_argument('-d', '--distance_cutoff', type=float, help=s, default=2.)
  args = parser.parse_args()
  assert os.path.exists(args.pdb_file)
  reduced_pdb_file = reduce_pdb(args.pdb_file)
  pairs = get_nonbonded_interactions(file_name=reduced_pdb_file,
                                     distance_cutoff=args.distance_cutoff)

  print 'Here are the inter-asymmetric-unit contacts :'
  for pair in pairs  :
    print pair

if __name__ == '__main__' :
  run(sys.argv[1:])

