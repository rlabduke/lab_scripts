#!/usr/bin/python

# Instructions from Liz:

# Run MolProbity on pdb id
#     - to generate xxxxH.pdb, multi-crit chart, and multi-crit kin with H-bonds, etc.

# Run cablam command in phenix:  $phenix.cablam_training cablam=True ~[path from phenix directory location to pdb_file_name]/3kwdH.pdb > 3kwdH_3CA_angle.txt
#      - to generate residue id info and the 3CA angle with residue CA in the middle (n-1,             n, n+1)

# Run DSSP on xxxx.pdb 
#     - to generate dssp.txt file (I have trimmed mine, but script could be revised to include                 'line.startswith')

# Run probe command:  $probe -once -mc -Radius[value ???  I used 1.0] -u -onedot "atom_HA_,atom_HA2 protein" "atom_O__ protein" xxxxH.pdb >> xxxxH_onedot_HACO.txt       **Note:  the two spaces following the O are correct.
#     - to generate text file including onedot (single entry for contacts) mcmc contacts for all                 HAs (including the HA2s on Gly) and C(O)s in a pdb file.

# Run Angle_add_HACO.py  (see attached)
#      - to combine pertinent data from all three previous runs into one text file and then transfer to an excel sheet for graphing, sorting.


echo "usage: HACO_data_onefile.bash filename 1>&2"

$muscle_pdb_path=/home/smlewis/whole_PDB_as_PDB/

#generate local copy of file
gunzip 


#Generate Re
