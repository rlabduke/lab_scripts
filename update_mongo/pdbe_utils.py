import sys
import json

urls = {
  'summary':"/pdb/entry/summary",
  'experiment':"/pdb/entry/experiment",
  'sifts':"/mappings"
}

PY3 = sys.version > '3'

if PY3:
    import urllib.request as urllib2
else:
    import urllib2

SERVER_URL = "https://www.ebi.ac.uk/pdbe/api"

def make_request(url, data):   
    request = urllib2.Request(url)

    try:
        url_file = urllib2.urlopen(request, data)
    except urllib2.HTTPError as e:
        if e.code == 404:
            print("[NOTFOUND %d] %s" % (e.code, url))
        else:
            print("[ERROR %d] %s" % (e.code, url))

        return None

    return url_file.read().decode()

def get_request(url, arg, pretty=False):
    full_url = "%s/%s/%s?pretty=%s" % (SERVER_URL, url, arg, str(pretty).lower())
    
    return make_request(full_url, None)

def post_request(url, data, pretty=False):
    full_url = "%s/%s/?pretty=%s" % (SERVER_URL, url, str(pretty).lower())
    
    if isinstance(data, (list, tuple)):
        data = ",".join(data)
    
    return make_request(full_url, data.encode())

class PDBdoc(object) :

  def __init__(self,pdb_id) :
    self.pdb_id = pdb_id
    self.document_types = urls.keys()

  def is_document_type(self,document_type) :
    return document_type in self.document_types

  def get_doc(self,document_type) :
    assert self.is_document_type(document_type)
    url = urls[document_type]
    pdbe_doc = get_request(url,self.pdb_id,pretty=True)
    if pdbe_doc == None : return
    pdbe_json_doc = json.loads(pdbe_doc)
    return pdbe_json_doc

  def write_doc(self,document_type,log=sys.stdout) :
    assert self.is_document_type(document_type)
    doc = self.get_doc(document_type)
    print >> log, json.dumps(doc,sort_keys=True,indent=4,
        separators=(',', ': '))

class Cluster(object) :

  def __init__(self,homology_level) :
    assert homology_level in [95,90,70,50]
    self.homology_level = homology_level

  def get_raw_data(self) :
    url = 'ftp://resources.rcsb.org/sequence/clusters/clusters%i.txt'
    data = str(self.homology_level)
    return make_request(url % self.homology_level,data)

  def get_cluster_dict(self) :
    rawdata = self.get_raw_data()
    if rawdata == None : return
    print type(rawdata)
    sl = rawdata.splitlines()
    print sl[:5]
