#!/usr/bin/python

import sys, os, os.path

usage = '''multiple_file_HACO.py list  Will run HACO_data_onefile.py on each PDB in list, putting results in own directories'''

if __name__ == '__main__':

    if (len(sys.argv) != 2) :
        print usage
        exit(1)

    #hardcoded for user videau on c3po
    HACO_data_onefile = "~/HACO_script/HACO_data_onefile.py "

    working_dir = os.getcwd()


    pdblist = sys.argv[1]
    for eachPDB in (open(pdblist, 'r').readlines()): #iterate through pdblist file
        eachPDB = eachPDB.rstrip() #strip trailing endline
        os.mkdir(eachPDB)
        os.chdir(eachPDB)
        HACO_cmd = HACO_data_onefile + eachPDB.upper()
        print HACO_cmd
        os.chdir(working_dir)
