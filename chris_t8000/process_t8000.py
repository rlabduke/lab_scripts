#! /usr/bin/python

import sys, os, os.path

usage = '''process_t8000.py $listfile          Will go look for PDBS from that listfile at (for 1UBQ) /home/smlewis/whole_PDB_as_PDB/ub/pdb1ubq.ent.gz, and grep out the chain listed into a new file'''


muscle_path = "/home/smlewis/whole_PDB_as_PDB/"

if __name__ == '__main__':

    if (len(sys.argv) != 2) :
        print usage
        exit(1)

    for eachLine in (open(sys.argv[1], 'r').readlines()):

        #extract the PDB name and chain id
        chainPDBPair = eachLine.split(',')
        pdb = chainPDBPair[0].lower() #NOTE the lower is an assumption for case-sensitive file systems
        chain = str(chainPDBPair[1]).rstrip() #clear we want string, for numerical chains

        #skip the header
        if (pdb == "pdb_id") :
            continue

        #path to this PDB
        interstice = pdb[1:3]
        inpath = muscle_path + interstice + "/pdb" + pdb + ".ent.gz"
        outpath = os.getcwd() + "/" +  pdb + "_" + chain + ".pdb"

        gunzip_subcommand = "gunzip --to-stdout " + inpath + " | "
        #this should grep the lines beginning with ATOM
        ATOM_subcommand = "egrep '^ATOM' | "
        #this will choose the lines with the right value in the chain column.  -F"" means column delimited.
        awk_subcommand = "awk -F '' '$22 == \"" + chain + "\"' "
        output_subcommand = "> " + outpath

        whole_command = gunzip_subcommand + ATOM_subcommand + awk_subcommand + output_subcommand
        print whole_command
        os.system(whole_command)
