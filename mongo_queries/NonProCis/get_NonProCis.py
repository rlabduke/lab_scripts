import os,sys
import argparse
sys.path.append('..')
import utils
import numpy

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
      elif isinstance(e,bool) : ls.append(str(e))
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
      # Skip if omegalyze doesn't exists.
      if not hasattr(mongores,'omegalyze') : continue
      for pres in mongores.prevres :
        if pres.resname not in utils.reslist : continue
        self.append(Omega(pres,mongores,self.resolution))
        #print type(mongores),type(pres)
      # Does it pass filters
      #break

  def write_csv(self,log=sys.stdout) :
    whead = True
    for i,omega in enumerate(self.omegas) :
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
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('-o','--homology_level', type=int,default=70,
                   help='Homology level can be 50, 70, 90, or 95. Default=70')
  parser.add_argument('-v','--verbose',action='store_true',help='Be verbose')
  args = parser.parse_args()

  # Get pdbs
  pdbs = utils.get_Top8000_pdb_list(homology_level=args.homology_level)
  if args.verbose :
    s = '%i pdbs found in Top8000 at homology level %i'
    utils.broadcast(s % (len(pdbs),args.homology_level))

  # Get connection to mongo
  mongocon = utils.MongodbConnection()

  # Iterate through pdbs
  counts = {"all_omega":
            {'n_unique':0,'n_alts':0,'n_unique_filter':0,'n_alts_filter':0},
          "pro":
            {'n_unique':0,'n_alts':0,'n_unique_filter':0,'n_alts_filter':0},
          "cispro":
            {'n_unique':0,'n_alts':0,'n_unique_filter':0,'n_alts_filter':0},
          "nonpro":
            {'n_unique':0,'n_alts':0,'n_unique_filter':0,'n_alts_filter':0},
          "nonprocis":
            {'n_unique':0,'n_alts':0,'n_unique_filter':0,'n_alts_filter':0}}
  nonprocis = NonProCis()
  for i,pc in enumerate(pdbs) :
    if i % 100 == 0 : print >> sys.stderr, "Through %i of %i.." % (i,len(pdbs))
    pdb_id,chain = pc
    #pdb_id,chain = "1avb","A"
    if args.verbose : print "working on %s %s..." % (pdb_id,chain)
    omegas = Omegas(pdb_id,chain,db=mongocon.db)
    #omegas.write_csv()
    omegas.set_counts()
    counts['all_omega']['n_unique'] += omegas.counts.all_omega_unique
    counts['all_omega']['n_unique_filter'] +=\
                              omegas.counts.all_omega_unique_filter
    counts['all_omega']['n_alts'] += omegas.counts.all_omega_alt
    counts['all_omega']['n_alts_filter'] += omegas.counts.all_omega_alt_filter
    counts['pro']['n_unique'] += omegas.counts.pro_unique
    counts['pro']['n_alts'] += omegas.counts.pro_alt
    counts['pro']['n_unique_filter'] += omegas.counts.pro_unique_filter
    counts['pro']['n_alts_filter'] += omegas.counts.pro_alt_filter
    counts['cispro']['n_unique'] += omegas.counts.cis_pro_unique
    counts['cispro']['n_alts'] += omegas.counts.cis_pro_alt
    counts['cispro']['n_unique_filter'] += omegas.counts.cis_pro_unique_filter
    counts['cispro']['n_alts_filter'] += omegas.counts.cis_pro_alt_filter
    counts['nonpro']['n_unique'] += omegas.counts.nonpro_unique
    counts['nonpro']['n_alts'] += omegas.counts.nonpro_alt
    counts['nonpro']['n_unique_filter'] += omegas.counts.nonpro_unique_filter
    counts['nonpro']['n_alts_filter'] += omegas.counts.nonpro_alt_filter
    counts['nonprocis']['n_unique'] += omegas.counts.cis_nonpro_unique
    counts['nonprocis']['n_alts'] += omegas.counts.cis_nonpro_alt
    counts['nonprocis']['n_unique_filter'] += omegas.counts.cis_nonpro_unique_filter
    counts['nonprocis']['n_alts_filter'] += omegas.counts.cis_nonpro_alt_filter
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

  print "In %i pdbs there were :" % i
  s = ': all omegas unique'
  v = '%i' % counts['all_omega']['n_unique']
  print v.ljust(10), s
  s = ': all omegas unique filter'
  v = '%.3f' % counts['all_omega']['n_unique_filter']
  print v.ljust(10), s
  s = ': all omegas alts'
  v = '%i' % counts['all_omega']['n_alts']
  print v.ljust(10), s
  s = ': all omegas alts filter'
  v = '%i' % counts['all_omega']['n_alts_filter']
  print v.ljust(10), s

  s = ': pro unique'
  v = '%i' % counts['pro']['n_unique']
  print v.ljust(10), s
  s = ': pro unique filter'
  v = '%.3f' % counts['pro']['n_unique_filter']
  print v.ljust(10), s
  s = ': pro alts'
  v = '%i' % counts['pro']['n_alts']
  print v.ljust(10), s
  s = ': pro alts filter'
  v = '%i' % counts['pro']['n_alts_filter']
  print v.ljust(10), s

  s = ': cispro unique'
  v = '%i' % counts['cispro']['n_unique']
  print v.ljust(10), s
  s = ': cispro unique filter'
  v = '%.3f' % counts['cispro']['n_unique_filter']
  print v.ljust(10), s
  s = ': cispro alts'
  v = '%i' % counts['cispro']['n_alts']
  print v.ljust(10), s
  s = ': cispro alts filter'
  v = '%i' % counts['cispro']['n_alts_filter']
  print v.ljust(10), s

  s = ': nonpro unique'
  v = '%i' % counts['nonpro']['n_unique']
  print v.ljust(10), s
  s = ': nonpro unique filter'
  v = '%.3f' % counts['nonpro']['n_unique_filter']
  print v.ljust(10), s
  s = ': nonpro alts'
  v = '%i' % counts['nonpro']['n_alts']
  print v.ljust(10), s
  s = ': nonpro alts filter'
  v = '%i' % counts['nonpro']['n_alts_filter']
  print v.ljust(10), s

  s = ': nonprocis unique'
  v = '%i' % counts['nonprocis']['n_unique']
  print v.ljust(10), s
  s = ': nonprocis unique filter'
  v = '%.3f' % counts['nonprocis']['n_unique_filter']
  print v.ljust(10), s
  s = ': nonprocis alts'
  v = '%i' % counts['nonprocis']['n_alts']
  print v.ljust(10), s
  s = ': nonprocis alts filter'
  v = '%i' % counts['nonprocis']['n_alts_filter']
  print v.ljust(10), s

if __name__ == '__main__' :
  run(sys.argv[1:])

