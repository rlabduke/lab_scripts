# HACO_script
script to automate Liz's HACO pipeline (carbon-hydrogen hbonds, usually in beta sheet)

There are three scripts, presented from innermost to outermost:

##Angle_add_HACO
Angle_add_HACO.py is a data multiplexer - it takes a CaBLAM file, a probe result, and a DSSP file and makes a nice data table.
Total commandline is "python Angle_add_HACO_02292016.py 3XXXH_cablam3CA_angle.txt 3XXX__dssp_1lineHeader.txt 3XXXH_onedot_1.0rad_All_02292016.txt PDBNAME"

##HACO_data_onefile
HACO_data_onefile.py takes a PDB id as an argument and generates the inputs needed for Angle_add_HACO (and then runs that sub script).  It is currently hardcoded for use as user videau on c3po, but with small changes can be set up for use on muscle instead.

It has labeled steps:
1. get the pdb
2. unzip
3. run reduce to strip hydrogens
4. run reduce to ADD hydrogens
5. delete the stripped pdb
6. run cablam on the reduced pdb
7. get dssp results from the web
8. run probe on the reduced pdb
9. run Angle_add_HACO

##multiple_file_HACO
multiple_file_HACO.py will take a list of PDB IDs (one per line) and run HACO_data_onefile on each, putting results in new directories.
