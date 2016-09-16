import unittest
import sys,os
assert os.getcwd().endswith('mongo_queries/tests')
sys.path.append('..')
import utils

class TestConnection(unittest.TestCase) :

  def __init__(self,*args,**kwargs) :
    super(TestConnection, self).__init__(*args, **kwargs)
    # Connect
    self.MDBConn = utils.MongodbConnection(db_name = 'pdb_info')

  def test_experiment(self) :
    cursor = self.MDBConn.db.experiment.find({"_id.pdb_id":"3KAT"})
    assert cursor.count() == 1
    self.assertEqual(cursor[0]['r_free_percent_reflections'],4.7)
    self.assertEqual(cursor[0]['percent_reflections_observed'],80.4)
    self.assertEqual(cursor[0]['r_factor'],0.202)
    self.assertEqual(cursor[0]['experimental_method'],'X-ray diffraction')

  def test_residues_colkeys(self) :
    q = {"pdb_id":"2efw","chain_id":"B","resseq":50}
    cursor = self.MDBConn.db.residues_colkeys.find(q)
    assert cursor.count() == 1
    rota = cursor[0]["rotalyze"]
    self.assertEqual(rota["rotamer"],"t80")
    self.assertEqual(rota["rotalyze_evaluation"],"Favored")
    self.assertEqual(rota["is_outlier"],False)
    self.assertEqual(cursor[0]["worst_all"]["twoFo_DFc_value"]["value"],
                     2.1679398668662495)

if __name__ == '__main__' :
    unittest.main()
