import os,sys
import argparse
sys.path.append('..')
import utils

class NQH_density(object) :

  def __init__(self,pdb_id) :
    self.pdb_id = pdb_id

  def set_density_values(self,db) :
    heads = ['pdb_id','chain_id','resseq','icode','altloc','resname','N','C_O']
    heads+= ["delta"]
    self.csv_lines = [heads]
    q = {"pdb_id":self.pdb_id,"resname":{"$in":["GLN","ASN","HIS"]}}
    cursor = db.residues_colkeys.find(q)
    for d in cursor :
      ll = [d["pdb_id"],d["chain_id"],"%i"%(d["resseq"]),d["icode"],d["altloc"]]
      ll+= [d["resname"]]
      if d["resname"] == "HIS" :
        Na = d["atoms"]["ND1"]["twoFo_DFc_value"]
        COa = d["atoms"]["CD2"]["twoFo_DFc_value"]
      elif d["resname"] == "GLN" :
        Na = d["atoms"]["NE2"]["twoFo_DFc_value"]
        COa = d["atoms"]["OE1"]["twoFo_DFc_value"]
      else :
        Na = d["atoms"]["ND2"]["twoFo_DFc_value"]
        COa = d["atoms"]["OD1"]["twoFo_DFc_value"]
      ll+= ["%.3f" % Na,"%.3f" % COa,"%.3f" % (Na-COa)]
      #print ",".join(ll)
      self.csv_lines.append(ll)

  def write_csv(self,out=sys.stdout) :
    for ll in self.csv_lines :
      print >> out, ",".join(ll)

def run(args) :
  desc = "Returns a csv of all map values ate atom positions of HQNs in the "
  desc+= "given PDB."
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdb_id', help='A pdb code')
  args = parser.parse_args()
  assert len(args.pdb_id) == 4

  mongocon = utils.MongodbConnection()
  nqhD = NQH_density(args.pdb_id)
  nqhD.set_density_values(mongocon.db)
  nqhD.write_csv()
# fastpig.r
if __name__ == '__main__' :
  run(sys.argv[1:])

