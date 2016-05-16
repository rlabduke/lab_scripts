import sys,os
from pymongo import MongoClient
import json

def get_creds_fn() :
  for dirname in sys.path:
    candidate = os.path.join(dirname, 'creds.json')
    if os.path.isfile(candidate):
      return candidate
  #raise RuntimeError("Can't find file creds.json")
  return None

def get_creds(user=None) :
  fn = get_creds_fn()
  if fn and os.path.exists(fn) :
    fle=open(fn,'r')
    d=json.load(fle)
    fle.close()
    return d['user'],d['dwp']
  import getpass
  if not user :
    user = getpass.getuser()
  print >> sys.stderr, "Please enter password for %s :" % user
  pwd = getpass.getpass()
  return user,pwd
collections = ['pdb_residues','file_info','experiment']

def print_json_pretty(d,log=sys.stdout) :
  try :
    print >> log, json.dumps(d,indent=4, separators=(',', ': '))
  except :
    print >> log, 'Could not print json document'

class MongodbConnection(object) :

  def __init__(self,db_name='pdb_info',user=None) :
    self.user, self.pwd = get_creds(user)
    self.db_name = db_name
    self.connect()

  def connect(self) :
    uri = "mongodb://%s:%s@daneel.research.duhs.duke.edu/"
    client = MongoClient(uri % (self.user,self.pwd))
    assert hasattr(client,self.db_name), '%s not on daneel'  % self.db_name
    self.db = getattr(client,self.db_name)