#!/usr/bin/python

import sys, os, os.path

usage = '''HACO_data_onefile.py PDBid          Will go look for the PDB with that id (for 1UBQ) /home/smlewis/whole_PDB_as_PDB/ub/pdb1ubq.ent.gz, and ????'''

muscle_path = "/home/smlewis/whole_PDB_as_PDB/"

if __name__ == '__main__':

    if (len(sys.argv) != 2) :
        print usage
        exit(1)

    pdb = sys.argv[1].lower() #NOTE the lower is an assumption for case-sensitive file systems
    cwd = os.getcwd() + "/"

    phenix_path = "/Applications/phenix-1.10.1-2155/build/setpaths.sh"        ##Added to source phenix automatically whenever this script is run, 
    load_phenix = "source " + phenix_path + " && "                            ##rather than using the 'do-phenix" command line alias command.
                                                                              ## "source" could also be included in the line18 string as "source /Applications/phenix-1.10.1-2155/build/setpaths.sh" 
    in_path = ""                                    ## Path to the alias do_phenix was found this way:
    muscle = False                                  ##  10:33:18 c3po ~> alias do_phenix
                                                    ##  source /Applications/phenix-1.10.1-2155/build/setpaths.csh
                                                    ##NOTE:  The alias path has ".csh" suffix, but the auto-source path is the same BUT has a ".sh" suffix.
    #path to this PDB
    if muscle:
        interstice = pdb[1:3]
        in_path = muscle_path + interstice + "/pdb" + pdb + ".ent.gz"
    else:
        in_path = "pdb" + pdb + ".ent.gz"
        #wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdbXXXX.ent.gz
        get_PDB_command = "curl ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb" + pdb + ".ent.gz > " + in_path
        print get_PDB_command
        os.system(get_PDB_command)


    #generate local copy of file
    unzipped_path = cwd + pdb + ".pdb"
    gunzip_command = "gunzip --to-stdout " + in_path + " > " + unzipped_path
    print gunzip_command
    os.system(gunzip_command)

    #run Reduce to strip old H  -  phenix is called, so the alias is auto-sourced with "load_phenix" as scripted above.
    reduce_stripped_path = cwd + pdb + ".stripped.pdb"
    reduce_strip_command = load_phenix + "phenix.reduce -quiet -trim -allalt " + unzipped_path + " | grep -v '^USER  MOD' > " + reduce_stripped_path
    print reduce_strip_command
    os.system(reduce_strip_command)

    #run Reduce  -  phenix is called, so the alias is auto-sourced with "load_phenix" as scripted above.
    reduced_path = cwd + pdb + "H.pdb"
    reduce_command = load_phenix + "phenix.reduce -quiet -nobuild9999 " + reduce_stripped_path + " > " + reduced_path
    print reduce_command
    os.system(reduce_command)

    #remove stripped pdb
    rm_strip_command = "rm " + reduce_stripped_path
    print rm_strip_command
    os.system(rm_strip_command)

    #run CaBLAM  -  phenix is called, so the alias is auto-sourced with "load_phenix" as scripted above.
    CaBLAM_results_path = cwd + pdb + "_3CA_angle.txt"
    CaBLAM_command = load_phenix + "phenix.cablam_training cablam=True " + reduced_path + " > " +  CaBLAM_results_path
    print CaBLAM_command
    os.system(CaBLAM_command)

    #get DSSP
    #this is dumb to use FTP, but when I tried rsyncing the whole database I realized they have literally every PDB's dssp in one directory, which causes performance issues
    #curl ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp/1ubq.dssp
    dssp_path = cwd + pdb + ".dssp"
    dssp_command = "curl ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp/" + pdb + ".dssp > " + dssp_path
    print dssp_command
    os.system(dssp_command)

    #cryptic remark
    # Run DSSP on xxxx.pdb
    #     - to generate dssp.txt file (I have trimmed mine, but script could be revised to include                 'line.startswith')

    #run Probe
    #**Note:  the two underscores following the O are correct.
    probe_path = cwd + pdb + "H_onedot_HACO.txt"
    ##probe_command = load_phenix + 'phenix.probe -once -mc -Radius1.0 -u -onedot "atom_HA_,atom_HA2 protein" "atom_O__ protein" ' + reduced_path + " > " + probe_path  ## added auto-source "load_phenix" in case this is ever needed.
    probe_command = '~/Desktop/Probe_2015/probe/probe -once -mc -Radius1.0 -u -onedot "atom_HA_,atom_HA2 protein" "atom_O__ protein" ' + reduced_path + " > " + probe_path
    print probe_command
    os.system(probe_command)    ## to change this script, pick one of the probe pathways, acccording to which probe version is most current?? best?? 

    # Run Angle_add_HACO.py to combine pertinent data from all three previous runs into one text file and then transfer to an excel sheet for graphing, sorting.
    #Total commandline is "python Angle_add_HACO_02292016.py 3XXXH_cablam3CA_angle.txt 3XXX__dssp_1lineHeader.txt 3XXXH_onedot_1.0rad_All_02292016.txt"
    HACO_path = cwd + pdb + ".HACOresult.csv"
    HACO_Liz_command = "python Angle_add_HACO.py " + CaBLAM_results_path + " " + dssp_path + " " + probe_path + " " + pdb + " > " + HACO_path
    print HACO_Liz_command
    os.system(HACO_Liz_command)
