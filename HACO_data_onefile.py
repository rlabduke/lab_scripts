#!/usr/bin/python

import sys, os, os.path

# Run DSSP on xxxx.pdb
#     - to generate dssp.txt file (I have trimmed mine, but script could be revised to include                 'line.startswith')

# Run probe command:  $probe -once -mc -Radius[value ???  I used 1.0] -u -onedot "atom_HA_,atom_HA2 protein" "atom_O__ protein" xxxxH.pdb >> xxxxH_onedot_HACO.txt       **Note:  the two spaces following the O are correct.
#     - to generate text file including onedot (single entry for contacts) mcmc contacts for all                 HAs (including the HA2s on Gly) and C(O)s in a pdb file.

# Run Angle_add_HACO.py  (see attached)
#      - to combine pertinent data from all three previous runs into one text file and then transfer to an excel sheet for graphing, sorting.

usage = '''HACO_data_onefile.py PDBid          Will go look for the PDB with that id (for 1UBQ) /home/smlewis/whole_PDB_as_PDB/ub/pdb1ubq.ent.gz, and ????'''

muscle_path = "/home/smlewis/whole_PDB_as_PDB/"

if __name__ == '__main__':

    if (len(sys.argv) != 2) :
        print usage
        exit(1)

    pdb = sys.argv[1].lower() #NOTE the lower is an assumption for case-sensitive file systems
    cwd = os.getcwd() + "/"

    #path to this PDB
    interstice = pdb[1:3]
    in_path = muscle_path + interstice + "/pdb" + pdb + ".ent.gz"
    unzipped_path = cwd + pdb + ".pdb"

    #generate local copy of file
    gunzip_command = "gunzip --to-stdout " + in_path + " > " + unzipped_path
    print gunzip_command

    #run Reduce
    reduced_path = cwd + pdb + "H.pdb"
    reduce_command = "phenix.reduce -quiet -trim -allalt " + unzipped_path + " | grep -v '^USER  MOD' > " + reduced_path
    print reduce_command

    #run CaBLAM
    CaBLAM_results_path = cwd + pdb + "_3CA_angle.txt"
    CaBLAM_command = "phenix.cablam_training cablam=True " + reduced_path + " > " +  CaBLAM_results_path
    print CaBLAM_command

    #get DSSP
    #this is dumb to use FTP, but when I tried rsyncing the whole database I realized they have literally every PDB's dssp in one directory, which causes performance issues
    #curl ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp/1ubq.dssp
    dssp_path = cwd + pdb + ".dssp"
    dssp_command = "curl ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp/" + pdb + ".dssp"
    print dssp_command


