import os,sys
import argparse
sys.path.append('..')
import utils
import numpy
import re
import pprint
import nqh_hingeCAdock as nqh_dock
from cctbx import geometry_restraints

dockatomlist = ["CA", "C", "O"]
ref_coords_dock = {"CA":[0,0,0],"C":[1.5,0,0],"O":[2.15,1,0]}

def translate(atomset, move):
  #print("move:"+str(move))
  #print("atomset:"+str(atomset))
  translated_atoms = {}
  for atom in atomset:
    if not atom is None and not atomset[atom] is None:
      translated_atoms[atom] = [ atomset[atom][0]-move[0], atomset[atom][1]-move[1], atomset[atom][2]-move[2] ]
  return translated_atoms

#move_coords and ref_coords are dicts with key atom_name -> value xyz coords in list form
def ca_dock(move_coords, resnum_to_dock="0"):
  #atoms_to_dock = get_atoms_to_flip(ag, resname, dock_atoms)
  #move_coords = {}
  #for atomname in dockatomlist:
  #  move_coords[atomname] = [all_atoms_to_dock[resnum_to_dock+atomname][0], all_atoms_to_dock[resnum_to_dock+atomname][1], all_atoms_to_dock[resnum_to_dock+atomname][2]]
  #print("move_coords, pre: ")
  #pprint.pprint(move_coords)
  #move_coords, pre: {'CA': [3.859, 46.816, 9.877], 'C': [4.85, 45.995, 9.061], 'O': [5.738, 46.544, 8.4]}
  translation_to_origin = move_coords[resnum_to_dock+dockatomlist[0]]
  move_coords = translate(move_coords, translation_to_origin)
  translation_to_origin = ref_coords_dock[dockatomlist[0]]
  #translate(ref_coords_dock, translation_to_origin)
  #print("move_coords, translated to origin: ")
  #pprint.pprint(move_coords)
  #return move_coords

  #CAdock calculation
  move_axis = [ move_coords[resnum_to_dock+dockatomlist[1]][0] - move_coords[resnum_to_dock+dockatomlist[0]][0],
                move_coords[resnum_to_dock+dockatomlist[1]][1] - move_coords[resnum_to_dock+dockatomlist[0]][1],
                move_coords[resnum_to_dock+dockatomlist[1]][2] - move_coords[resnum_to_dock+dockatomlist[0]][2] ]
  ref_axis =  [  ref_coords_dock[dockatomlist[2]][0] -  ref_coords_dock[dockatomlist[0]][0],
                 ref_coords_dock[dockatomlist[2]][1] -  ref_coords_dock[dockatomlist[0]][1],
                 ref_coords_dock[dockatomlist[2]][2] -  ref_coords_dock[dockatomlist[0]][2] ]
  #print("move_axis:"+str(move_axis))
  #print("ref_axis:"+str(ref_axis))
  
  alignaxis = numpy.cross(move_axis, ref_axis)
  #may need safe handling of cmag here
  #  probably not due to guarantee of distance >>0
  dotp = numpy.dot(move_axis,ref_axis)
  normalize1 = numpy.sqrt( move_axis[0]**2 + move_axis[1]**2 + move_axis[2]**2 )
  normalize2 = numpy.sqrt( ref_axis[0]**2 + ref_axis[1]**2 + ref_axis[2]**2 )
  angle1 = numpy.arccos( dotp/(normalize1*normalize2) ) #in radians, of course
  
  for atomname in move_coords:
    move_coords[atomname] = nqh_dock.doaxisrot(move_coords[atomname], angle1, alignaxis)
  #for atomname in all_atoms_to_dock:
  #  all_atoms_to_dock[atomname] = nqh_dock.doaxisrot(move_coords[atomname], angle1, alignaxis)
  #First rotation of CA dock complete
  #print("move_coords, after first rotation: ")
  #pprint.pprint(move_coords)
  #return move_coords
  #start second rotation of CA dock
  
  #get a dihedral
  #m3,r1,r2,r3
  angle2 = geometry_restraints.dihedral(sites=[move_coords[resnum_to_dock+dockatomlist[2]], ref_coords_dock[dockatomlist[0]], ref_coords_dock[dockatomlist[2]], ref_coords_dock[dockatomlist[1]]], angle_ideal=-40, weight=1).angle_model
  #in degrees, convert to radians
  angle2 = numpy.deg2rad(angle2)
  ##print >> sys.stderr, angle2
  
  for atomname in move_coords:
    move_coords[atomname] = nqh_dock.doaxisrot(move_coords[atomname], angle2, ref_axis)
  #for atomname in all_atoms_to_dock:
  #  all_atoms_to_dock[atomname] = nqh_dock.doaxisrot(move_coords[atomname], angle1, alignaxis)


  #translate from origin to original position
  #translate(move_coords, [-translation_to_origin[0], -translation_to_origin[1], -translation_to_origin[2]])

  #CA dock complete
  #print("move_coords, post: ")
  #pprint.pprint(move_coords)
  return move_coords

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

  model_num = 1
  file_num = 0
  fragment_set = {}
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
      if mongores.is_outlier() or not mongores.passes_filter('bb'):
        if args.verbose:
          print("excluding "+str(mongores)+" because it has outlier")
      else:
        # prevres and nextres are lists because those residues might have alts
        if mongores.altloc=="":
          if len(mongores.prevres) == 1:
            prev_residue = mongores.prevres[0]
            if prev_residue.is_outlier() or not prev_residue.passes_filter('bb') or not prev_residue.altloc=="":
              if args.verbose:
                print("excluding "+str(mongores)+" because previous residue has outlier or alt")
            else:
              if len(mongores.nextres) == 1:
                next_residue = mongores.nextres[0]
                if next_residue.is_outlier() or not next_residue.passes_filter('bb') or not next_residue.altloc=="":
                  if args.verbose:
                    print("excluding "+str(mongores)+" because next residue has outlier or alt")
                else:
                  #sys.stdout.write('.')
                  #sys.stdout.flush()
                  #print("building fragment from "+str(mongores))
                  if model_num == 9999:
                    model_num = 1
                    file_num = file_num + 1
                    out_file.close()
                    print(fragment_set)
                    #assert False
                    out_file = open("output_fragments"+str(file_num)+".pdb", 'w')
                  fragment = utils.MongoPdbFragment([prev_residue, mongores, next_residue])
                  fragment.set_bb_atoms(ca_dock(fragment.get_bb_atoms()))
                  if len(fragment_set) == 0: 
                    fragment_set[fragment] = 1
                  else:
                    rmsd = 1
                    for test_frag in fragment_set:
                      rmsd = fragment.get_rmsd(test_frag)
                      if rmsd <= 0.5:
                        break
                    if rmsd > 0.5:
                      fragment_set[fragment] = 1
                      #print(fragment.get_bb_atoms())
                      #print(fragment.get_atom_records('bb'))
                      #assert False
                      out_file.write("MODEL{:>9}\n".format(model_num))
                      out_file.write(fragment.get_atom_records(translated=True, region='bb'))
                      out_file.write("ENDMDL\n")
                      model_num = model_num+1
                      #print("res "+str(mongores)+" prev res: "+str(prev_residue)+" next res: "+str(next_residue))
                    else:
                      fragment_set[test_frag] = fragment_set[test_frag] + 1
  out_file.close()

if __name__ == '__main__' :
  run(sys.argv[1:])

