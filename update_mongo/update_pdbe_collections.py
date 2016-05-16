import os,sys
import argparse
import pdbe_utils
import mongo_utils
import json

collections = ['experiment','summary']

def get_current_pdbs() :
  pdblistfn = 'allpdbs.l'
  assert os.path.exists(pdblistfn)
  pdbs = {}
  n = 0
  fle = open(pdblistfn,'r')
  for l in fle :
    if not l.strip().endswith('.ent.gz') : continue
    pdb_id = l[l.rfind('/pdb')+4:l.rfind('.ent.gz')].upper()
    assert len(pdb_id) == 4
    mid = pdb_id[1:3]
    if mid not in pdbs.keys() : pdbs[mid] = []
    pdbs[mid].append(pdb_id)
    n += 1
  return n,pdbs

def get_missing(pdbs,pdbs_in_db) :
  missing = []
  for mid,lst in pdbs.items() :
    for pdb_id in lst :
      if mid not in pdbs_in_db.keys() or pdb_id not in pdbs_in_db[mid] :
         missing.append(pdb_id)
      #if pdb_id not in pdbs_in_db[mid] : missing.append(pdb_id)
      #if pdb_id in pdbs_in_db[mid] : print pdb_id,"in"
      #	else : print pdb_id,"out" 
  return missing

def get_existing_pdb_ids(collection,mcon) :
  assert collection in collections
  assert hasattr(mcon.db,collection)
  dbcol = getattr(mcon.db,collection)
  cursor = dbcol.find({},{"_id":True})
  pdb_ids = {}
  i = 0
  for i,_id in enumerate(cursor) :
    #print _id
    #if collection == 'experiment' : pdb_id = _id['_id']["pdb_id"]
    #else : pdb_id = _id['_id']
    pdb_id = _id['_id']["pdb_id"]
    mid = pdb_id[1:3]
    if mid not in pdb_ids : pdb_ids[mid] = []
    if pdb_id not in pdb_ids[mid] : pdb_ids[mid].append(pdb_id)
    #if i > 5 : break
  if i == 0 : ln = i
  else : ln = i + 1
  return ln,pdb_ids

def run(args) :
  desc = """Updates the either pdb_info.experiment or pdb_info.summary on 
  daneel. The script has the following steps :

  - get a list of currently deposited pdbs (in allpdbs.l) 
  - get a list of pdbs in given collection. 
  - get a list of deposited pdbs missing from given collection
  - iterate trough the missing pdbs :
    - get the given pdb's experiment documenti[s] from the PDBe
    - insert that/those document into given collection 

"""
  parser = argparse.ArgumentParser(description=desc,
                                 formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('-v','--verbose',help='Verbose outpur',
                      action='store_true')
  parser.add_argument('-t','--collection', type=str,
                      help='experiment or summary',default='experiment')
  args = parser.parse_args()
  assert args.collection in ['experiment','summary']


  # get a list of currently deposited pdbs (in allpdbs.l) 
  n_pdbs,pdbs = get_current_pdbs()
  print >> sys.stderr, 'There are currently %i deposited pdbs.' % n_pdbs

  # connect to mongo on daneel
  mcon = mongo_utils.MongodbConnection()

  # get a list of pdbs in given collection
  n_pdbs_in_db,pdbs_in_db = get_existing_pdb_ids(args.collection,mcon)
  s = 'There are currently %i pdbs in the collection "%s".'
  print >> sys.stderr, s % (n_pdbs_in_db,args.collection)

  # get a list of deposited pdbs missing from pdb_info.experiment
  missing = get_missing(pdbs,pdbs_in_db)
  s = 'There are currently %i pdbs missing from in the collection "%s".'
  print >> sys.stderr, s % (len(missing),args.collection)
  #exit()

  # iterate trough the missing pdbs
  if len(missing) > 3000 : factor = 1000
  else : factor = 100
  print >> sys.stderr, '\n\nBegin iterating missing pdbs...\n'
  msg = '\n%i records inserted - %.2f %% done.\n'
  for i,pdb_id in enumerate(missing) :
    if args.verbose : print >> sys.stderr, "working on %s..." % pdb_id
    # get the given pdb's document from the PDBe
    pdoc = pdbe_utils.PDBdoc(pdb_id)
    doc = pdoc.get_doc(args.collection)
    #print doc
    #exit()
    if doc is None :
      s = '\nWARNING: Skipping %s -- PDBe request returned None\n' % pdb_id
      print >> sys.stderr, s
      continue
    # I am assuming that the documents are a dict of length 1 and the key
    # is the pdb id with the value being a list of one or more dict (which 
    # is/are the one[s] we're after). But I am unsure that this is universal 
    # thus the foloowing.
    #print_json_pretty(doc)
    assert len(doc) == 1, doc
    assert len(doc.keys()) == 1,len(doc.keys())
    k = doc.keys()[0]
    assert k.upper() == pdb_id
    assert type(doc[k]) is list,type(doc[k])
    assert type(doc[k][0]) is dict,type(doc[k][0])
    #print_json_pretty(doc[k][0])

    # insert that/those document into pdb_info.experiment
    for j,mdoc in enumerate(doc[k]) :
      mdoc["_id"] = {"pdb_id":pdb_id,"n":j}
      dbcol = getattr(mcon.db,args.collection)
      dbcol.insert(mdoc)
      if args.verbose :
        print >> sys.stderr, '  %s inserted into pdb_info.%s' % (pdb_id,args.collection)
      if i%factor == 0 : print >> sys.stderr, msg % (i,(i*100.0)/len(missing))
    if i > 5 : break
  print >> sys.stderr, '\n\n%i total records inserted\n\n' % i

if __name__ == '__main__' :
  run(sys.argv[1:])

