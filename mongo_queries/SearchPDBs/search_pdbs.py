import os,sys
import argparse
scriptpath = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(scriptpath,'..'))
import utils

default_resol = 1.2
exmethods = ["X-ray diffraction","X-ray solution scattering",
  "X-ray powder diffraction","Solution NMR","Solid-state NMR",
  "Electron Microscopy","Neutron Diffraction","Infrared spectroscopy",
  "Electron crystallography","Fiber diffraction","Fluorescence transfer"]

resolution_methods = ["X-ray diffraction","Electron Microscopy",
  "Neutron Diffraction","Electron crystallography","Fiber diffraction"]

def get_summary_query(args) :
  q = {}
  if args.dna : q["number_of_entities.dna"] = {"$gt":0}
  if args.rna : q["number_of_entities.rna"] = {"$gt":0}
  if args.protein : q["number_of_entities.protein"] = {"$gt":0}
  if args.xdna : q["number_of_entities.dna"] = 0
  if args.xrna : q["number_of_entities.rna"] = 0
  if args.xprotein : q["number_of_entities.protein"] = 0
  return q

def get_experiment_query(args) :
  q = {}
  if args.experimental_method != '0' :
    q["experimental_method"] = args.experimental_method
  if q["experimental_method"] in resolution_methods :
    if args.resolution > 0 : q["resolution"] = {"$lte":args.resolution}
  if not args.agnostic_structure_factors :
    q["experiment_data_available"] = "Y"
  if args.spacegroup :
    s = "spacegroup can only be spcified when experimental_method is"
    s+= " X-ray diffraction."
    assert q["experimental_method"] == "X-ray diffraction",s
    q["spacegroup"] = args.spacegroup
  return q

def get_pdbs(cursor) :
  pdbs = []
  for d in cursor :
    #print d
    pdb_id = str(d["_id"]["pdb_id"])
    assert len(pdb_id) == 4
    if not pdb_id in pdbs : pdbs.append(pdb_id)
    #break
  #exit() 
  return pdbs

def get_experiment_pdbs(args,mongocon) :
  query = get_experiment_query(args)
  utils.print_json_pretty(query,log=sys.stderr)
  cursor = mongocon.db.experiment.find(query,{"_id":1})
  return get_pdbs(cursor)

def get_summary_pdbs(args,mongocon) :
  query = get_summary_query(args)
  if query == {} : return
  utils.print_json_pretty(query,log=sys.stderr)
  cursor = mongocon.db.summary.find(query,{"_id":1})
  return get_pdbs(cursor)

def run() :
  # Parse arguments
  desc = """This script is meant to return lists of pdbs that fit the given criteria. The criteria are set using flags. To see all flags set the help flage, -h. Seeing as you are already reading this, you likely already set the gelp flag. Good job!"""
  parser = argparse.ArgumentParser(description=desc)
  # experiment parameters
  s = 'Get pdbs with resolution higher than this. Default is %.1f. To turm off'
  s+= ' this filter set to 0 or less .'
  parser.add_argument('-e','--resolution',type=float,default=default_resol,
                      help= s % default_resol)
  parser.add_argument('-s','--spacegroup',type=str,
                      help="Return PDBs with this spacegroup")
  parser.add_argument('-f','--agnostic_structure_factors',action='store_true',
           help='Return PDBs regardless of availability of structure factors')
  # summary parameters
  parser.add_argument('-d','--dna',action='store_true',
                      help='Get structures with dna')
  parser.add_argument('-r','--rna',action='store_true',
                      help='Get structures with rna')
  parser.add_argument('-p','--protein',action='store_true',
                      help='Get structures with protein')
  parser.add_argument('--xdna',action='store_true',
                      help='Exclude structures with dna')
  parser.add_argument('--xrna',action='store_true',
                      help='Exclude structures with rna')
  parser.add_argument('--xprotein',action='store_true',
                      help='Exclude structures with protein')
  s = 'Get structures with this experimental method. Default is "X-ray '
  s+= 'diffraction". To turn off this filter set to 0. Possible choices are '
  s+= '"0", "%s"' % '", "'.join(exmethods)
  parser.add_argument('-m','--experimental_method', type=str,
                      default="X-ray diffraction", help=s)
  args = parser.parse_args()

  # We have to query two databases, experiment and summary, to get all the
  # required parameters.

  # connect to mongo
  mongocon = utils.MongodbConnection()

  # get list of PDBs that match the given experment parameters
  experiment_pdbs = get_experiment_pdbs(args,mongocon)
  #print len(experiment_pdbs)
  #print experiment_pdbs[:8]

  # get list of PDBs that match the given summary parameters
  summary_pdbs = get_summary_pdbs(args,mongocon)
  #print len(summary_pdbs)
  #print summary_pdbs[:8]

  # get the union
  if summary_pdbs :
    pdbs = []
    for pdb_id in summary_pdbs :
      if pdb_id in experiment_pdbs : pdbs.append(pdb_id)
  else : pdbs = experiment_pdbs

  log = sys.stdout
  for pdb_id in pdbs :
    print >> log, pdb_id

  s = "%i PDBs fetched"
  print >> sys.stderr, s % len(pdbs)

if __name__ == "__main__" :
  run()
