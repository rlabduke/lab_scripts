import os,sys
import argparse

def get_pdbs_with_sf() :
  print >> sys.stderr, "Getting pdbs with structure factors..."
  sf_pdbs = []
  count = 0
  rootd = '/net/cci/pdb_mirror/structure_factors/'
  for (dirpath, dirnames, filenames) in os.walk(rootd):
    for dr in dirnames :
      for pat2,dirs2,fls2 in os.walk(os.path.join(dirpath,dr)) :
        for fn in fls2 :
          assert fn.startswith('r'), os.path.join(pat2,fn)
          pdb_id = fn[1:5]
          sf_pdbs.append(fn[1:5])
          count += 1
          if count % 10000 == 0 : print 'through %i pdbs' % count
  print >> sys.stderr, "There are %i pdbs with structure factors\n" % count
  return sf_pdbs

def get_pdbs_in_db(pdbs_in_db_fn) :
  print >> sys.stderr, "Getting pdbs in db"
  pdbs_in_db = {}
  fle = open(pdbs_in_db_fn,'r')
  for i,line in enumerate(fle) :
    pdb = line.strip()
    mid = pdb[1:3]
    if mid not in pdbs_in_db.keys() : pdbs_in_db[mid] = []
    pdbs_in_db[mid].append(pdb)
    if i % 10000 == 0 : print >> sys.stderr, 'through %i pdbs' % i
  fle.close()
  print >> sys.stderr, "There are %i pdbs in the db\n" % i
  return pdbs_in_db

def get_pdbs_to_be_done(pdbs_in_db,pdbs_with_sf,out) :
  print >> sys.stderr, "Getting pdbs to be done..."
  count = 0
  for i,pdb in enumerate(pdbs_with_sf) :
    mid = pdb[1:3]
    if pdb not in pdbs_in_db[mid] :
      print >> out, pdb
      count += 1
    if i % 10000 == 0 : print >> sys.stderr, 'through %i pdbs' % i
  print >> sys.stderr, "There are %i pdbs that need to be done\n" % count

def run(args) :
  desc = "This script creates a pdbs_to_be_done.l file needed to run "
  desc+= "run_mongo_validation.py on the LBL cluster. You must be on the LBL "
  desc+= "cluster to run this script as it uses the PDB mirror at LBL to get "
  desc+= "the current PDB codes (with structure factors)." 
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdbs_in_db_fn', help='A list of PDBs already in the DB.')
  args = parser.parse_args()
  assert os.path.exists(args.pdbs_in_db_fn)

  # get PDBs with structure factors in the PDB mirror
  pdbs_with_sf = get_pdbs_with_sf()

  # get pdbs in db 
  pdbs_in_db = get_pdbs_in_db(args.pdbs_in_db_fn)
  # write pdbs_to_be_done.l
  outfile = 'pdbs_to_be_done.l'
  fle = open(outfile,'w')
  get_pdbs_to_be_done(pdbs_in_db   = pdbs_in_db,
                      pdbs_with_sf = pdbs_with_sf,
                      out          = fle)
  fle.close()
  print >> sys.stderr, "\n%s written.\n" % outfile

if __name__ == '__main__' :
  run(sys.argv[1:])
