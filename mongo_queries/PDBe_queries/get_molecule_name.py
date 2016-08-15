import os,sys
import argparse
sys.path.append('..')
from utils import pdbe_utils

class Entity(object) :

  def __init__(self,pdb_id,entityd) :
    assert type(entityd) is dict
    self.pdb_id = pdb_id
    self.entityd = entityd
    self.molecule_type = entityd["molecule_type"]

  def get_info(self) :
    return {'pdb_id':self.pdb_id,
            'molecule_type':self.entityd["molecule_type"],
	    'entity_id':self.entityd["entity_id"],
            'in_chains':','.join(self.entityd["in_chains"]),
            'molecule_name':','.join(self.entityd["molecule_name"])}

  def write_info(self,write_head=False,log=sys.stdout) :
    info = self.get_info()
    heads = ['pdb_id','entity_id','in_chains','molecule_type','molecule_name']
    if write_head : print >> log, ':'.join(heads)
    print >> log, '\t'.join([str(info[label]) for label in heads])

def get_molecules(pdb_id) :
  doc = pdbe_utils.PDBdoc(pdb_id)
  jsond = doc.get_doc("molecule")
  molecules = []
  for entity in jsond[pdb_id] :
    molecules.append(Entity(pdb_id,entity))
  return molecules

def write_molecule_info_csv(pdb_id,
                            molecule_type='polypeptide',
                            log=sys.stdout) :
  molecules = get_molecules(pdb_id)
  assert type(molecules) is list
  assert len(molecules) > 0
  for molecule in molecules :
    mt = molecule.molecule_type
    if mt=='all' or mt.startswith(molecule_type) : 
      molecule.write_info()

def run(args) :
  desc = "This script takes either a PDB code or a list of codes and returns "
  desc+= "the molecule name[s] of said PDB to stdout. The list of PDBs is "
  desc+= "expected to be a file where each line is a PDB code and just a PDB "
  desc+= "code."
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdblist_or_code',
                      help='A file of PDB codes or just a PDB code.',
                      type=str)
  args = parser.parse_args()
  pdbs_done = []
  if os.path.exists(args.pdblist_or_code) :
    fle = open(args.pdblist_or_code)
    for i,l in enumerate(fle) :
      #if i > 10 : break
      if i%100 == 0 : print >> sys.stderr, "Trough %i" % i
      pdb_id = l.strip()
      if pdb_id in pdbs_done : continue
      pdbs_done.append(pdb_id)
      write_molecule_info_csv(pdb_id)
    fle.close()
  else :
    pdb_id = args.pdblist_or_code
    assert len(pdb_id) == 4
    write_molecule_info_csv(pdb_id)

if __name__ == '__main__' :
  run(sys.argv[1:])

