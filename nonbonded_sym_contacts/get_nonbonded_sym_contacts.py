import os,sys
import argparse
import libtbx.load_env
from mmtbx import monomer_library
from mmtbx.monomer_library import pdb_interpretation
from libtbx import group_args
from libtbx import easy_run
import re
import iotbx.pdb

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
    self.set_residue_ids()

  def __str__(self) :
    return "%s - %.2f - %s"% \
        (self.residues[0].id,self.distance,self.residues[1].id)

  def set_pair_labels(self,i_seq,j_seq) :
    ilab = parse_label(self.scatterers[i_seq].label)
    jlab = parse_label(self.scatterers[j_seq].label)
    self.residues = (ilab,jlab)

  def get_residue_id(reslf,res) :
    l = [res.alt_loc,res.resname,res.chain,res.resseq,res.icode]
    return ":".join(l)

  def set_residue_ids(self) :
    r1 = self.get_residue_id(self.residues[0])
    r2 = self.get_residue_id(self.residues[1])
    self.residue_ids = (r1,r2)

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

class NonbondedIinteractions(list) :

  def get_nonbonded_interactions(self,
                                 file_name,
                                 distance_cutoff = 2, 
                                 select = "nucleotide",
                                 filter_residues = (None,None),
                                 filter_atoms = (None,None)) :
    print >> sys.stderr, "Source residues : %s" % filter_residues[0]
    print >> sys.stderr, "Target residues : %s" % filter_residues[1]
    assert select in ["nucleotide",None]
    # filter_residues must be a 2tuple where each element must be None or list
    assert isinstance(filter_residues,tuple)
    assert len(filter_residues) == 2
    if filter_residues[0] : assert isinstance(filter_residues[0],list)
    if filter_residues[1] : assert isinstance(filter_residues[1],list)
    self.file_name = file_name
    # filter_atoms must be a 2tuple where each element must be None or dict, if
    # dict the keys must match the residues in filter_residues
    assert isinstance(filter_atoms,tuple)
    assert len(filter_atoms) == 2
    if filter_atoms[0] :
      assert isinstance(filter_atoms[0],dict)
      assert filter_residues[0]
      for k in filter_residues[0] : assert k in filter_atoms[0],k
    if filter_atoms[1] :
      assert isinstance(filter_atoms[1],dict)
      assert filter_residues[1]
      for k in filter_residues[1] : assert k in filter_atoms[1],k
    # get geostd
    chem_data = libtbx.env.find_in_repositories(
      relative_path="chem_data/geostd",
      test=os.path.isdir)
    # setup monomer_library server
    if chem_data is not None:
      mon_lib_srv = monomer_library.server.server()
      ener_lib = monomer_library.server.ener_lib()
    # get NA selection
    pdb_inp = iotbx.pdb.input(file_name=self.file_name)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    if select :
      asc = pdb_hierarchy.atom_selection_cache()
      selection = asc.selection("nucleotide")
      pdb_hierarchy = pdb_hierarchy.select(selection)
      pdb_hierarchy.atoms().reset_i_seq()
    # processed pdb file
    pdb_processed_file = pdb_interpretation.process(
      #file_name=self.file_name,
      pdb_inp=pdb_hierarchy.as_pdb_input(),
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      crystal_symmetry=pdb_inp.crystal_symmetry(),
      force_symmetry = True)
    #print dir(pdb_processed_file)
    #exit()
    # get grm
    grm = pdb_processed_file.geometry_restraints_manager(
      hard_minimum_nonbonded_distance=0,
      show_energies = False,
      plain_pairs_radius = 5.0)
    xrs = pdb_processed_file.xray_structure()
    #xrs = xrs_all#.select(selection)
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
    while i < n_nonb and sorted_nonb[i][3] < distance_cutoff :
      # skip intra-asym interactions
      if sorted_nonb[i][-1] is None : i += 1; continue
      pair = nonbondedPair(sorted_nonb[i],xrs.scatterers())
      rn0 = pair.residues[0].resname
      rn1 = pair.residues[1].resname
      an0 = pair.residues[0].name
      an1 = pair.residues[1].name
      if filter_residues[0] and rn0 not in filter_residues[0] : i += 1; continue
      if filter_residues[1] and rn1 not in filter_residues[1] : i += 1; continue
      if filter_atoms[0] and an0 not in filter_atoms[0][rn0] : i += 1; continue
      if filter_atoms[1] and an1 not in filter_atoms[1][rn1] : i += 1; continue
      self.append(pair)
      i += 1

  def write_raw_pairs(self,log=sys.stdout) :
    for pair in self : print >> log, pair

  def write_formatted_pairs(self,log=sys.stdout) :
    found = '%i pair[s] found in %s'
    residue_level_pairs = {}
    found_n = 0
    for pair in self :
      rn1 = pair.residues[0].resname
      rn2 = pair.residues[1].resname
      #print pair
      key = (pair.residue_ids[0],pair.residue_ids[1])
      if not key in residue_level_pairs.keys() : residue_level_pairs[key] = []
      t = (pair.residues[0].name,pair.residues[1].name,'%.2f'%pair.distance)
      residue_level_pairs[key].append(t)
    if len(residue_level_pairs) > 0 :
      print >> sys.stderr, found % (len(residue_level_pairs),self.file_name)
      found_n += len(residue_level_pairs)
    for k,l in residue_level_pairs.items() :
      print >> log, k
      for e in l : print >> log, '    ' + str(e)
    #break
    print >> sys.stderr, '%i TOTAL pair[s] found' % found_n

def parseargs() :
  default_distance = 2.
  desc = "This script will show you inter-asymmetric-unit contacts from a PDB "
  desc+= "with a CRYST card. I suggest 3g6t as it has a nice base pair between "
  desc+= " asymmetric units, C 1 to C 2. This script adds hysrogens for you. "
  desc+= "\n\nAGAIN, THESE ARE ONLY INTER_ASYMMETRIC_UNIT INTERACTIONS!!\n\n"
  parser = argparse.ArgumentParser(description=desc,
                                  formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('pdb_file', help='A pdb file')
  s = 'Distance cutoff between the two interacting atoms. default is %.1f'
  parser.add_argument('-d', '--distance_cutoff', type=float,
                       help=s % default_distance, default=default_distance)
  s = 'Write a residue-formatted output'
  parser.add_argument('-f','--formatted',action="store_true",help=s)
  return parser.parse_args()

def run() :
  args = parseargs()
  assert os.path.exists(args.pdb_file)
  reduced_pdb_file = reduce_pdb(args.pdb_file)
  pairs = NonbondedIinteractions()
  pairs.get_nonbonded_interactions(file_name=reduced_pdb_file,
                                   distance_cutoff=args.distance_cutoff,
                                   select = None,
                                   filter_residues = (None,None))

  s = 'Here are the inter-asymmetric-unit contacts within a distance of %.1f :'
  print s % args.distance_cutoff
  if args.formatted : pairs.write_formatted_pairs()
  else : pairs.write_raw_pairs()
  # cleanup
  os.remove(reduced_pdb_file)

if __name__ == '__main__' :
  run()

