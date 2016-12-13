import os,sys
import argparse
sys.path.append('..')
import utils

class NQH_density(object) :

  def __init__(self,pdb_id) :
    self.pdb_id = pdb_id

  def set_density_values(self,db) :
    #heads = ['pdb_id','chain_id','resseq','icode','altloc','resname']
    #heads+= ['worst_rscc','worst_2fo-fc','N','C','O','dN','dC','dO','delta']
    self.csv_lines = []
    q = {"pdb_id":self.pdb_id,"resname":{"$in":["GLN","ASN","HIS"]}}
    cursor = db.residues_colkeys.find(q)
    for d in cursor :
      ll = [d["pdb_id"],d["chain_id"],"%i"%(d["resseq"]),d["icode"].strip()]
      ll+= [d["altloc"],d["resname"],"%.3f" % d["worst_all"]["rscc"]["value"]]
      ll+= ["%.3f" % d["worst_all"]["twoFo_DFc_value"]["value"]]
      if d["resname"] == "HIS" :
        if not "ND1" in d["atoms"].keys() and not "CD2" in d["atoms"].keys() :
          continue
        Na = d["atoms"]["ND1"]["twoFo_DFc_value"] or ''
        Ca = d["atoms"]["CD2"]["twoFo_DFc_value"]
        Oa = ''
        dNa = d["atoms"]["ND1"]["Fo_DFc_value"] or ''
        dCa = d["atoms"]["CD2"]["Fo_DFc_value"]
        dOa = ''
      elif d["resname"] == "GLN" :
        if not "NE2" in d["atoms"].keys() and not "OE1" in d["atoms"].keys() :
          continue
        Na = d["atoms"]["NE2"]["twoFo_DFc_value"]
        Ca = ''
        Oa = d["atoms"]["OE1"]["twoFo_DFc_value"]
        dNa = d["atoms"]["NE2"]["Fo_DFc_value"]
        dCa = ''
        dOa = d["atoms"]["OE1"]["Fo_DFc_value"]
      else :
        if not "ND2" in d["atoms"].keys() and not "OD1" in d["atoms"].keys() :
          continue
        Na = d["atoms"]["ND2"]["twoFo_DFc_value"]
        Ca = ''
        Oa = d["atoms"]["OD1"]["twoFo_DFc_value"]
        dNa = d["atoms"]["ND2"]["Fo_DFc_value"]
        dCa = ''
        dOa = d["atoms"]["OD1"]["Fo_DFc_value"]
      if d["resname"] == "HIS" :
        ll+= ["%.3f" % Na,"%.3f" % Ca,Oa]
        ll += ["%.3f" % (Na-Ca)]
        ll+= ["%.3f" % dNa,"%.3f" % dCa,dOa]
      else :
        ll+= ["%.3f" % Na,Ca,"%.3f" % Oa]
        ll += ["%.3f" % (Oa-Na)]
        ll+= ["%.3f" % dNa,dCa,"%.3f" % dOa]
      #print ",".join(ll)
      self.csv_lines.append(ll)

  def write_csv(self,out=sys.stdout) :
    heads = ['pdb_id','chain_id','resseq','icode','altloc','resname']
    heads+= ['worst_rscc','worst_2fo-fc','N','C','O','dN','dC','dO','delta']
    exp = "# delta is the 2fo-fc values for the O - N for GLN and ASN and "
    exp+= "N - C for HIS."
    print >> out, exp
    print >> out, ",".join(heads)
    self.csv_lines.sort(key=lambda x: int(x[2]))
    for ll in self.csv_lines :
      assert len(heads) == len(ll)
      print >> out, ",".join(ll)

def run(args) :
  desc = "Returns a csv of all map values at atom positions of HQNs in the "
  desc+= "given PDB."
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdb_id', help='A pdb code')
  args = parser.parse_args()
  assert len(args.pdb_id) == 4

  mongocon = utils.MongodbConnection()
  nqhD = NQH_density(args.pdb_id)
  nqhD.set_density_values(mongocon.db)
  nqhD.write_csv()

if __name__ == '__main__' :
  run(sys.argv[1:])

