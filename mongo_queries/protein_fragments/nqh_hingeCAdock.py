from __future__ import division
import os, sys
from iotbx import pdb
from cctbx import geometry_restraints
import numpy
import math

ROTATOMLIST = {
  "ASN":[' CB ',' CG ',' OD1',' ND2','','HD21','HD22'],
  "GLN":[' CG ',' CD ',' OE1',' NE2','HE21','HE22'],
  "HIS":[' CB ',' CG ',' CE1',' NE2',' ND1',' CD2',' HD1',' HE1',' HD2',' HE2']}

HINGEATOMLIST = {
  "ASN":[' CG ',' OD1',' ND2','','HD21','HD22'],
  "GLN":[' CD ',' OE1',' NE2','HE21','HE22'],
  "HIS":[' CG ',' CE1',' NE2',' ND1',' CD2',' HD1',' HE1',' HD2',' HE2']}

DOCKATOMLIST = {
  "ASN":[' CA ',' OD1',' ND2',' CB ',' CG ','','HD21','HD22',' HB2',' HB3'],
  "GLN":[' CA ',' OE1',' NE2',' CB ',' CG ',' CD ','HE21','HE22',' HB2',' HB3',' HG2',' HG3'],
  "HIS":[' CA ',' CE1',' NE2',' CB ',' ND1',' CD2',' CG ',' HD1',' HE1',' HD2',' HE2',' HB2',' HB3']}

def vectorize(p1, p2):
  v = [ p2[0]-p1[0] , p2[1]-p1[1] , p2[2]-p1[2] ]
  return v

#Returns the scalar length of a vector
def veclen(v):
  return math.sqrt( v[0]**2 + v[1]**2 + v[2]**2 )

def parse_user_mods(filename):
  flipped_residues = []
  f = file(filename, 'rb')
  for line in f.readlines():
    if line.startswith("USER  MOD"):
      if "FLIP" in line:
        temp = line.strip().split(':')
        flipped_residues.append(temp[1])
  f.close()
  return flipped_residues

def get_atoms_to_flip(ag,resname,needed_atoms):
  fetched_atoms= {}
  atomlist = needed_atoms[resname]
  rg = ag.parent()
  alts = [ag.altloc, ' ', '']
  for atom in rg.atoms():
    if atom.name in atomlist and atom.id_str()[9:10] in alts:
      fetched_atoms[atom.name] = atom
  return fetched_atoms

def translate(atomset, move):
  for atom in atomset:
    atomset[atom] = [ atomset[atom][0]-move[0], atomset[atom][1]-move[1], atomset[atom][2]-move[2] ]

#rotates move_atom around alignaxis by angle
#This is the real mathematical workhorse, pretty much everythibng else is selecting atoms or determining alignaxis and angle
#I think this is "just" as rotation matix
#             [x,y,z]?, radians, [x,y,z]?
def doaxisrot(move_atom, angle, alignaxis, a=[0.0,0.0,0.0]):
  xxyyzz = [ alignaxis[0] - a[0], alignaxis[1] - a[1], alignaxis[2] - a[2] ]
  div = numpy.sqrt(xxyyzz[0]**2 + xxyyzz[1]**2 + xxyyzz[2]**2)
  cosn = [ xxyyzz[0]/div, xxyyzz[1]/div, xxyyzz[2]/div ]

  costheta = numpy.cos(angle)
  sintheta = numpy.sin(angle)

  a11 = cosn[0]*cosn[0] + (1-cosn[0]*cosn[0])*costheta
  a12 = cosn[0]*cosn[1]*(1-costheta) + cosn[2]*sintheta
  a13 = cosn[0]*cosn[2]*(1-costheta) - cosn[1]*sintheta

  a21 = cosn[0]*cosn[1]*(1-costheta) - cosn[2]*sintheta
  a22 = cosn[1]*cosn[1] + (1-cosn[1]*cosn[1])*costheta
  a23 = cosn[1]*cosn[2]*(1-costheta) + cosn[0]*sintheta

  a31 = cosn[0]*cosn[2]*(1-costheta) + cosn[1]*sintheta
  a32 = cosn[1]*cosn[2]*(1-costheta) - cosn[0]*sintheta
  a33 = cosn[2]*cosn[2] + (1-cosn[2]*cosn[2])*costheta

  fx1 = move_atom[0]-alignaxis[0]
  fy1 = move_atom[1]-alignaxis[1]
  fz1 = move_atom[2]-alignaxis[2]

  fx2 = fx1*a11 + fy1*a21 + fz1*a31
  fy2 = fx1*a12 + fy1*a22 + fz1*a32
  fz2 = fx1*a13 + fy1*a23 + fz1*a33

  return [ fx2+alignaxis[0] , fy2+alignaxis[1] , fz2+alignaxis[2] ]

#Takes a selected sidechain (ag in hierarchy terms) and performs a series of transformations on it
#First, a 180 degree rotation (one rotation)
#Second, a hinge motion in case the terminal group was non-planar (one rotation)
#Third, a three-point dock using the CA and the two terminal atoms of the sidechain (two rotation)
def do_rot_then_flip(ag, rot_atoms=ROTATOMLIST, hinge_atoms=HINGEATOMLIST, dock_atoms=DOCKATOMLIST):
  resname = ag.id_str()[1:4]
  dockatomlist = dock_atoms[resname]
  hingeatomlist = hinge_atoms[resname]
  rotatomlist = rot_atoms[resname]
  atoms_to_rot = get_atoms_to_flip(ag, resname, ROTATOMLIST)
  atoms_to_hinge = get_atoms_to_flip(ag, resname, HINGEATOMLIST)
  atoms_to_dock = get_atoms_to_flip(ag, resname, DOCKATOMLIST)
  move_coords, ref_coords_hinge, ref_coords_dock = {}, {}, {}
  for atomname in atoms_to_rot:
    move_coords[atomname] = [atoms_to_rot[atomname].xyz[0], atoms_to_rot[atomname].xyz[1], atoms_to_rot[atomname].xyz[2]]
  #reference coordinates for the original position are establiched now, since each step will update the model
  for atomname in atoms_to_hinge:
    ref_coords_hinge[atomname] = [atoms_to_hinge[atomname].xyz[0], atoms_to_hinge[atomname].xyz[1], atoms_to_hinge[atomname].xyz[2]]
  for atomname in atoms_to_dock:
    ref_coords_dock[atomname] = [atoms_to_dock[atomname].xyz[0], atoms_to_dock[atomname].xyz[1], atoms_to_dock[atomname].xyz[2]]

  #start 180 degree rotation
  #Each step starts with a translation to the origin
  translation_to_origin = move_coords[rotatomlist[0]]
  translate(move_coords, translation_to_origin)
  #Align axis is the CB-CG vector (CG-CD for GLN)
  alignaxis = [ move_coords[rotatomlist[1]][0] - move_coords[rotatomlist[0]][0],
                move_coords[rotatomlist[1]][1] - move_coords[rotatomlist[0]][1],
                move_coords[rotatomlist[1]][2] - move_coords[rotatomlist[0]][2] ]

  for atomname in move_coords:
    move_coords[atomname] = doaxisrot(move_coords[atomname], numpy.pi, alignaxis)

  #translate from origin to original position
  translate(move_coords, [-translation_to_origin[0], -translation_to_origin[1], -translation_to_origin[2]])
  change_coords(atoms_to_rot, move_coords)
  #180 degree rotation complete

  
  #hinge motion setup
  atoms_to_hinge = get_atoms_to_flip(ag, resname, hinge_atoms)
  move_coords = {}
  for atomname in atoms_to_hinge:
    move_coords[atomname] = [atoms_to_hinge[atomname].xyz[0], atoms_to_hinge[atomname].xyz[1], atoms_to_hinge[atomname].xyz[2]]

  translation_to_origin = ref_coords_hinge[hingeatomlist[0]]
  translate(move_coords, translation_to_origin)
  translate(ref_coords_hinge, translation_to_origin)

  #hinge motion calculation
  #lots of cross products, basically
  #Calculate normal to the plane of the original sidechain and normal to the rotated sidechain
  normplaneold = numpy.cross(vectorize(ref_coords_hinge[hingeatomlist[0]],ref_coords_hinge[hingeatomlist[1]]),vectorize(ref_coords_hinge[hingeatomlist[1]],ref_coords_hinge[hingeatomlist[2]]))
  normplanenew = numpy.cross(vectorize(move_coords[hingeatomlist[1]],move_coords[hingeatomlist[2]]),vectorize(move_coords[hingeatomlist[0]],move_coords[hingeatomlist[1]]))
  #axis of rotation is the line of intersection between the planes of the old and new sidechains
  #due to translation to the origin, this is given by the cross product of the plane normals
  hingeaxis = numpy.cross(normplaneold,normplanenew)
  #angle of rotation alge between the planes of the old and new sidechains
  #this is found using the dot product of the plane normals
  hingeangle = -1*numpy.arccos(numpy.dot(normplaneold,normplanenew) / (veclen(normplaneold)*veclen(normplanenew)))

  for atomname in move_coords:
    move_coords[atomname] = doaxisrot(move_coords[atomname], hingeangle, hingeaxis)
  
  #translate from origin to original position
  translate(move_coords, [-translation_to_origin[0], -translation_to_origin[1], -translation_to_origin[2]])
  change_coords(atoms_to_hinge, move_coords)
  #hinge motion complete

  #CAdock setup
  atoms_to_dock = get_atoms_to_flip(ag, resname, dock_atoms)
  move_coords = {}
  for atomname in atoms_to_dock:
    move_coords[atomname] = [atoms_to_dock[atomname].xyz[0], atoms_to_dock[atomname].xyz[1], atoms_to_dock[atomname].xyz[2]]

  #move_coords and ref_coords are dicts with key atom_name -> value xyz coords in list form
  print("move_coords:"+str(move_coords))
  print("ref_coords_doc:"+str(ref_coords_dock))
  print("docatomlist:"+str(dockatomlist))
  translation_to_origin = ref_coords_dock[dockatomlist[0]]
  translate(move_coords, translation_to_origin)
  translate(ref_coords_dock, translation_to_origin)
  print("move_coords, translate to origin: "+str(move_coords))

  #CAdock calculation
  move_axis = [ move_coords[dockatomlist[1]][0] - move_coords[dockatomlist[0]][0],
                move_coords[dockatomlist[1]][1] - move_coords[dockatomlist[0]][1],
                move_coords[dockatomlist[1]][2] - move_coords[dockatomlist[0]][2] ]
  ref_axis =  [  ref_coords_dock[dockatomlist[2]][0] -  ref_coords_dock[dockatomlist[0]][0],
                 ref_coords_dock[dockatomlist[2]][1] -  ref_coords_dock[dockatomlist[0]][1],
                 ref_coords_dock[dockatomlist[2]][2] -  ref_coords_dock[dockatomlist[0]][2] ]
  print("move_axis:"+str(move_axis))
  print("ref_axis:"+str(ref_axis))

  alignaxis = numpy.cross(move_axis, ref_axis)
  #may need safe handling of cmag here
  #  probably not due to guarantee of distance >>0
  dotp = numpy.dot(move_axis,ref_axis)
  normalize1 = numpy.sqrt( move_axis[0]**2 + move_axis[1]**2 + move_axis[2]**2 )
  normalize2 = numpy.sqrt( ref_axis[0]**2 + ref_axis[1]**2 + ref_axis[2]**2 )
  angle1 = numpy.arccos( dotp/(normalize1*normalize2) ) #in radians, of course
  
  for atomname in move_coords:
    move_coords[atomname] = doaxisrot(move_coords[atomname], angle1, alignaxis)
  #First rotation of CA dock complete
  
  #start second rotation of CA dock
  
  #get a dihedral
  #m3,r1,r2,r3
  angle2 = geometry_restraints.dihedral(sites=[move_coords[dockatomlist[2]], ref_coords_dock[dockatomlist[0]], ref_coords_dock[dockatomlist[2]], ref_coords_dock[dockatomlist[1]]], angle_ideal=-40, weight=1).angle_model
  #in degrees, convert to radians
  angle2 = numpy.deg2rad(angle2)
  ##print >> sys.stderr, angle2
  
  for atomname in move_coords:
    move_coords[atomname] = doaxisrot(move_coords[atomname], angle2, ref_axis)

  #translate from origin to original position
  translate(move_coords, [-translation_to_origin[0], -translation_to_origin[1], -translation_to_origin[2]])

  #update hierarchy with new coordinates
  change_coords(atoms_to_dock, move_coords)
  #CA dock complete


def change_coords(atoms_to_flip, move_coords):
  for atomname in atoms_to_flip:
    #This is the hierarchy way to change coordinates
    atoms_to_flip[atomname].set_xyz(move_coords[atomname])


if __name__ == '__main__':
  file_to_flip = sys.argv[1]
  flip_records = sys.argv[2]
  
  flipped_residues = parse_user_mods(flip_records)
  
  pdb_io = pdb.input(file_to_flip)
  hierarchy = pdb_io.construct_hierarchy()
  
  flip_rg_ids = []
  flip_ag_ids = []
  
  for flip in flipped_residues:
    chain_id = flip[0:2]
    resseq = flip[2:6]
    ins = flip[6:7]
    resname = flip[7:10]
    altloc = flip[14:15]
    flip_rg_ids.append(chain_id+resseq+ins)
    flip_ag_ids.append(altloc + resname +chain_id + resseq + ins)
  
  sourcefile = open(file_to_flip)
  for line in sourcefile:
    #print line[0:5]
    if line.startswith("CRYST"):
      sys.stdout.write(line)
  sourcefile.close()
  
  #rg.id_str() returns ccnnnni
  #ag.id_str() returns arrrccnnnni
  for ag in hierarchy.atom_groups():
    if ag.id_str() in flip_ag_ids:
      do_rot_then_flip(ag)
#temp removing this for debugging
#  print hierarchy.as_pdb_string()
  
  #the stderr printing is largely for my ability to track which residues have been flipped.
  for flip in flipped_residues:
    printid = flip[0:2]+flip[2:6]+flip[6:7]+flip[14:15]
    #       =      cc  +     nnnn+     i   +     a
    print >> sys.stderr, printid#flip
  #for rg in hierarchy.residue_groups():
  #  if rg.id_str() in flip_rg_ids:
    

