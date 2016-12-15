import os,sys
import argparse
sys.path.append('..')
import utils
import numpy
import re

class NonProCis(list) :

  def add_nonprocis(self) : pass

class Omega(object) :

  def __init__(self,mr0,mr1,resolution) :
    assert isinstance(mr0,utils.MongoResidue)
    assert isinstance(mr1,utils.MongoResidue)
    self.residue0 = mr0
    self.residue1 = mr1
    self.resname = mr1.resname
    self.resolution = resolution
    self.residue0_passes_filter = self.residue0.passes_filter(region='bb')
    self.residue1_passes_filter = self.residue1.passes_filter(region='bb')
    self.omega = self.residue1.omegalyze.omega
    self.omega_type = self.residue1.omegalyze.type

  def get_phipsi(self,num) :
    assert num in [0,1], num
    res = getattr(self,'residue%i' % num)
    if not hasattr(res,'ramalyze') : return None,None
    return res.ramalyze.phi, res.ramalyze.psi

  def passes_filter(self) :
    return self.residue0_passes_filter and self.residue1_passes_filter

  def is_non_pro_cis(self) :
    if self.resname != 'PRO' and self.omega_type == 'Cis' : return True
    else : return False

  def as_csv(self,withheads=False,log=sys.stdout) :
    heads = ['res0','res0_pass_filter','res1','res1_pass_filter']
    heads+= ['resolution','omega','omega_type','phi0','psi0','phi1','psi1']
    phi0,phsi0 = self.get_phipsi(0)
    phi1,phsi1 = self.get_phipsi(1)
    if withheads : print >> log, ','.join(heads)
    l = [self.residue0.as_str(),self.residue0_passes_filter]
    l+= [self.residue1.as_str(),self.residue1_passes_filter]
    l+= [self.resolution,self.omega,self.omega_type]
    l+= [phi0,phsi0,phi1,phsi1]
    ls = []
    for e in l :
      if e == None : ls.append("None")
      elif isinstance(e,float) : ls.append('%.2f' % e)
      elif isinstance(e,bool) or isinstance(e,int) : ls.append(str(e))
      else : ls.append(e)
    print >> log, ','.join(ls)

class Omegas(object) :

  def __init__(self,pdb_id,chain,db) :
    self.pdb_id = pdb_id
    self.chain = chain
    self.db = db
    self.omegas = []
    self.set_omegas()

  def __len__(self) : return len(self.omegas)

  def append(self,omega) :
    assert isinstance(omega,Omega)
    self.omegas.append(omega)

  def set_omegas(self) :
    residues = utils.MongoResidueList(db     = self.db,
                                      pdb_id = self.pdb_id,
                                      chain  = self.chain)
    self.resolution = residues.get_resolution()
    # iterate through residues
    keys = residues.ordered_keys()
    for k in keys :
      mongores = residues[k]
      if mongores.resname not in utils.reslist : continue
      # Skip : return False if omegalyze doesn't exists.
      if not hasattr(mongores,'omegalyze') : continue
      for pres in mongores.prevres :
        if pres.resname not in utils.reslist : continue
        self.append(Omega(pres,mongores,self.resolution))
        #print type(mongores),type(pres)
      # Does it pass filters
      #break

  def write_csv(self,noncis_pro_only=True,log=sys.stdout) :
    whead = False
    for i,omega in enumerate(self.omegas) :
      if noncis_pro_only and not omega.is_non_pro_cis() : continue
      omega.as_csv(withheads=whead,log=log)
      if i == 0 : whead = False

  def set_counts(self) :
    self.counts = utils.group_args(all_omega_unique = 0,##
                                   all_omega_alt = 0,##
                                   all_omega_unique_filter = 0,##
                                   all_omega_alt_filter = 0,##
                                   pro_unique = 0,#
                                   pro_alt = 0,#
                                   pro_unique_filter = 0,#
                                   pro_alt_filter = 0,#
                                   cis_pro_unique = 0,#
                                   cis_pro_alt = 0,#
                                   cis_pro_unique_filter = 0,#
                                   cis_pro_alt_filter = 0,#
                                   nonpro_unique = 0,##
                                   nonpro_alt = 0,##
                                   nonpro_unique_filter = 0,#
                                   nonpro_alt_filter = 0,#
                                   cis_nonpro_unique = 0,#
                                   cis_nonpro_alt = 0,#
                                   cis_nonpro_unique_filter = 0,#
                                   cis_nonpro_alt_filter = 0,#
                                   )
    uniqes = []
    uniqes_filter = []
    for omega in self.omegas :
      noalt = omega.residue1.as_str(noalt=True)
      passes = omega.passes_filter()
      occ = omega.residue1.lowest_occ()
      # ispro
      if omega.residue1.resname == 'PRO' : ispro = True
      else : ispro = False
      # iscis
      if omega.residue1.omegalyze.type == "Cis" : iscis = True
      else : iscis = False
      # start couting
      self.counts.all_omega_alt += 1
      if passes :
        self.counts.all_omega_alt_filter += 1
        self.counts.all_omega_unique_filter += 1*occ
      if ispro :
        self.counts.pro_alt += 1
        if passes :
          self.counts.pro_alt_filter += 1
          self.counts.pro_unique_filter += 1*occ
      else :
        self.counts.nonpro_alt += 1
        if passes :
          self.counts.nonpro_alt_filter += 1
          self.counts.nonpro_unique_filter += 1*occ
      # cis stuff
      if iscis :
        if ispro :
          self.counts.cis_pro_alt += 1
          if passes :
            self.counts.cis_pro_alt_filter += 1
            self.counts.cis_pro_unique_filter += 1*occ
        else :
          self.counts.cis_nonpro_alt += 1
          if passes :
            self.counts.cis_nonpro_alt_filter += 1
            self.counts.cis_nonpro_unique_filter += 1*occ
        
      # is it uniqe?
      if noalt not in uniqes :
        uniqes.append(noalt)
        self.counts.all_omega_unique += 1
        if ispro :
          self.counts.pro_unique += 1
          if iscis : self.counts.cis_pro_unique += 1
        else :
          self.counts.nonpro_unique += 1
          if iscis : self.counts.cis_nonpro_unique += 1

def has_non_pro_cis(pdb_id,chain,db) :
  q = {"pdb_id":pdb_id,"chain":chain,"resname":{"$ne":"PRO"}}
  q["omegalyze.type"] = "CIS"
  if db.residues_colkeys.count(q) == 0 : return False
  return True

def get_cis_residues(residues) :
  residue_keys = residues.ordered_keys()
  npc_residues  = []
  for reskey in residue_keys :
    mongores = residues[reskey]
    if not hasattr(mongores,'omegalyze') : continue
    if mongores.omegalyze.type == "Cis" :
      npc_residues.append(mongores)
  return npc_residues

def get_npc_residues(residues) :
  residue_keys = residues.ordered_keys()
  npc_residues  = []
  for reskey in residue_keys :
    mongores = residues[reskey]
    if not hasattr(mongores,'omegalyze') : continue
    if mongores.omegalyze.type == "Cis" and mongores.resname != 'PRO' :
      npc_residues.append(mongores)
  return npc_residues

def run(args) :
  desc = "A query script to getnon-pro cis residues from the Top8000 at a "
  desc+= "given homology level."
  print(desc)
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('-o','--homology_level', type=int,default=70,
                   help='Homology level can be 50, 70, 90, or 95. Default=70')
  parser.add_argument('-v','--verbose',action='store_true',help='Be verbose')
  args = parser.parse_args()

  # Get connection to mongo
  mongocon = utils.MongodbConnection()

  # Get pdbs
  pdbs = utils.get_Top8000_pdb_list(homology_level=args.homology_level,
                                    verbose=True,
                                    connection=mongocon)
  if args.verbose :
    s = '%i pdbs found in Top8000 at homology level %i'
    utils.broadcast(s % (len(pdbs),args.homology_level))

  # Set noncispro filename
  #ncp_fn = 'non_cis_pro_filtered_%i.csv' % args.homology_level
  #ncp_log = open(ncp_fn,'w')
  #heads = ['res0','res0_pass_filter','res1','res1_pass_filter']
  #heads+= ['resolution','omega','omega_type','phi0','psi0','phi1','psi1']
  #print >> ncp_log, ','.join(heads)

  # Iterate through pdbs
  #counts = {"all_omega":
  #          {'n_unique':0,'n_alts':0,'n_unique_filter':0,'n_alts_filter':0},
  #        "pro":
  #          {'n_unique':0,'n_alts':0,'n_unique_filter':0,'n_alts_filter':0},
  #        "cispro":
  #          {'n_unique':0,'n_alts':0,'n_unique_filter':0,'n_alts_filter':0},
  #        "nonpro":
  #          {'n_unique':0,'n_alts':0,'n_unique_filter':0,'n_alts_filter':0},
  #        "nonprocis":
  #          {'n_unique':0,'n_alts':0,'n_unique_filter':0,'n_alts_filter':0}}
  #nonprocis = NonProCis()
  model_num = 1
  file_num = 0
  out_file = open("output_fragments0.pdb", 'w')
  mongocon.set_db(db='pdb_info')
  for i,pc in enumerate(pdbs) :
    if i % 100 == 0 : print >> sys.stderr, "Through %i of %i.." % (i,len(pdbs))
    pdb_id,chain = pc
    #pdb_id,chain = "1d3g","A"
    if args.verbose : print "working on %s %s..." % (pdb_id,chain)
    
    #print(mongocon.db)
    residues = utils.MongoResidueList(mongocon.db, pdb_id, chain)
    #print(residues)
    keys = residues.ordered_keys()
    #print(keys)
    for k in keys :
      mongores = residues[k]
      if mongores.is_outlier() and mongores.passes_filter('bb'):
        if args.verbose:
          print("excluding "+str(mongores)+" because it has outlier")
      else:
        # prevres and nextres are lists because those residues might have alts
        if mongores.altloc=="":
          if len(mongores.prevres) == 1:
            prev_residue = mongores.prevres[0]
            if prev_residue.is_outlier() and mongores.passes_filter('bb'):
              if args.verbose:
                print("excluding "+str(mongores)+" because previous residue has outlier")
            else:
              if len(mongores.nextres) == 1:
                next_residue = mongores.nextres[0]
                if next_residue.is_outlier() and mongores.passes_filter('bb'):
                  if args.verbose:
                    print("excluding "+str(mongores)+" because next residue has outlier")
                else:
                  if model_num == 5001:
                    model_num = 1
                    file_num = file_num + 1
                    out_file.close()
                    out_file = open("output_fragments"+str(file_num)+".pdb", 'w')
                  out_file.write("MODEL{:>9}\n".format(model_num))
                  out_file.write(prev_residue.get_atom_records('bb'))
                  out_file.write(mongores.get_atom_records('bb'))
                  out_file.write(next_residue.get_atom_records('bb'))
                  out_file.write("ENDMDL\n")
                  model_num = model_num+1
                  #print("res "+str(mongores)+" prev res: "+str(prev_residue)+" next res: "+str(next_residue))
  out_file.close()
#print(keys);
#_id,model_id,altloc,chain_id,atoms,pdb_id,resseq,icode,restype,resname,worst_all,worst_sc,omegalyze,ramalyze,worst_bb,cablam,rotalyze
      
    #omegas = Omegas(pdb_id,chain,db=mongocon.db)
    #omegas.set_counts()
    #counts['all_omega']['n_unique'] += omegas.counts.all_omega_unique
    #counts['all_omega']['n_unique_filter'] +=\
    #                          omegas.counts.all_omega_unique_filter
    #counts['all_omega']['n_alts'] += omegas.counts.all_omega_alt
    #counts['all_omega']['n_alts_filter'] += omegas.counts.all_omega_alt_filter
    #counts['pro']['n_unique'] += omegas.counts.pro_unique
    #counts['pro']['n_alts'] += omegas.counts.pro_alt
    #counts['pro']['n_unique_filter'] += omegas.counts.pro_unique_filter
    #counts['pro']['n_alts_filter'] += omegas.counts.pro_alt_filter
    #counts['cispro']['n_unique'] += omegas.counts.cis_pro_unique
    #counts['cispro']['n_alts'] += omegas.counts.cis_pro_alt
    #counts['cispro']['n_unique_filter'] += omegas.counts.cis_pro_unique_filter
    #counts['cispro']['n_alts_filter'] += omegas.counts.cis_pro_alt_filter
    #counts['nonpro']['n_unique'] += omegas.counts.nonpro_unique
    #counts['nonpro']['n_alts'] += omegas.counts.nonpro_alt
    #counts['nonpro']['n_unique_filter'] += omegas.counts.nonpro_unique_filter
    #counts['nonpro']['n_alts_filter'] += omegas.counts.nonpro_alt_filter
    #counts['nonprocis']['n_unique'] += omegas.counts.cis_nonpro_unique
    #counts['nonprocis']['n_alts'] += omegas.counts.cis_nonpro_alt
    #counts['nonprocis']['n_unique_filter'] += omegas.counts.cis_nonpro_unique_filter
    #counts['nonprocis']['n_alts_filter'] += omegas.counts.cis_nonpro_alt_filter
    # write non-cis pro csv if they exist
    #if omegas.counts.cis_nonpro_unique_filter > 0 : omegas.write_csv(log=ncp_log)
    # get all_aa counts
    #counts['all_aa']['n_unique'] += residues.counts.unique_canonical_aa
    #counts['all_aa']['n_unique_filter'] += \
    #                 residues.counts.unique_canonical_aa_filter
    #counts['all_aa']['n_alts'] += residues.counts.all_canonical_aa
    #counts['all_aa']['n_alts_filter'] += residues.counts.all_canonical_aa_filter
    # get non-pro cis residues
    #cis_residues = get_cis_residues(residues)
    # Sanity check
    #assert len(cis_residues) > 0
    #break
    #if i > 15 : break

  #ncp_log.close()
  #print >> sys.stderr, '%s written.' % ncp_fn
  #re_fn = 'counts_%i.csv' % args.homology_level
  #re_log = open(re_fn,'w')
  #print >> re_log, "In %i pdbs there were :" % i
  #s = ': all omegas unique'
  #v = '%i' % counts['all_omega']['n_unique']
  #print >> re_log, v.ljust(10), s
  #s = ': all omegas unique filter'
  #v = '%.3f' % counts['all_omega']['n_unique_filter']
  #print >> re_log, v.ljust(10), s
  #s = ': all omegas alts'
  #v = '%i' % counts['all_omega']['n_alts']
  #print >> re_log, v.ljust(10), s
  #s = ': all omegas alts filter'
  #v = '%i' % counts['all_omega']['n_alts_filter']
  #print >> re_log, v.ljust(10), s

  #s = ': pro unique'
  #v = '%i' % counts['pro']['n_unique']
  #print >> re_log, v.ljust(10), s
  #s = ': pro unique filter'
  #v = '%.3f' % counts['pro']['n_unique_filter']
  #print >> re_log, v.ljust(10), s
  #s = ': pro alts'
  #v = '%i' % counts['pro']['n_alts']
  #print >> re_log, v.ljust(10), s
  #s = ': pro alts filter'
  #v = '%i' % counts['pro']['n_alts_filter']
  #print >> re_log, v.ljust(10), s

if __name__ == '__main__' :
  run(sys.argv[1:])

