import os, sys
from libtbx import easy_run

outd = 'val_out'
pgm_path = '/net/marbles/raid1/bhintze/PDB_mongodb/run_mongo_validation.py'

def get_pdb_code(pat) :
  assert os.path.exists(pat)
  return os.path.basename(pat)[:4]

def get_pdb_paths() :
  fn = 'pdbs_to_be_done.l'
  paths = []
  pa = '/net/cci/pdb_mirror/mmcif/%s/%s.cif.gz'
  fle = open(fn,'r')
  for l in fle :
    ls = l.strip().lower()
    if ls == '' : continue
    assert len(ls) == 4
    mid = ls[1:3]
    pat = pa % (mid,ls)
    assert os.path.exists(pat), pat
    paths.append(pat)
  fle.close()
  return paths

def get_cmd(mmcif_fp) :
  #pdbid = get_pdb_code(mmcif_fp)
  pdbid = mmcif_fp
  cmd_ = ["phenix.python",pgm_path,pdbid,'--outdir',outd,'--pdb_cif',mmcif_fp]
  #cmd_ = ["phenix.python",pgm_path,pdbid,'--outdir',outd]
  cmd_ += ['--write_out_file','-d','residue']
  return ' '.join(cmd_)

assert os.path.exists(outd)
pdb_paths = get_pdb_paths()

#for i, pdb_path in enumerate(pdb_paths):
#  pdb_code = get_pdb_code(pdb_path)
#  cmd = get_cmd(pdb_path)
#  print cmd
#  if i > 4 : break
#sys.exit()
#print get_pdb_code(pdb_paths[0])
#print get_cmd(pdb_paths[0])
#pdb_paths = pdb_paths[:9]
#print len(pdb_paths);sys.exit()

def run(only_i=None,
        chunk_n=1,
        chunk_size=len(pdb_paths)
        ):
  try: only_i = int(only_i)
  except ValueError: only_i=None

  try: chunk_n = int(chunk_n)
  except ValueError: chunk_n=1
  try: chunk_size = int(chunk_size)
  except ValueError: chunk_size=len(pdb_paths)

  assert chunk_size==len(pdb_paths) or chunk_n==1
  assert chunk_n>0
  assert chunk_size>0
  if chunk_n!=1:
    chunk_size = (len(pdb_paths)-1)//chunk_n+1
  elif chunk_size!=1:
    chunk_n = len(pdb_paths)%chunk_size+1

  for i, pdb_path in enumerate(pdb_paths):
    pdb_code = get_pdb_code(pdb_path)
    #pdb_code = pdb_path
    cmd = get_cmd(pdb_path)
    if only_i is None:
      print i, cmd
      continue
    else:
      if chunk_n!=1 or chunk_size!=len(pdb_paths):
        if only_i!=i//chunk_size+1: continue
      else:
        if only_i!=i+1: continue
    print 'Running', i+1
    print 'Command', cmd
    easy_run.call(cmd)

if __name__=="__main__":
  run(*tuple(sys.argv[1:]))
