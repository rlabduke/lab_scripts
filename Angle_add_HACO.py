#!/usr/bin/python

import os, sys
#Total commandline is "python Angle_add_HACO_02292016.py 3XXXH_cablam3CA_angle.txt 3XXX__dssp_1lineHeader.txt 3XXXH_onedot_1.0rad_All_02292016.txt"
# sys.argv is ["python Angle_add_HACO_02292016.py", "cablam_file", "dssp_file", "Probe_file"]
# sys.argv[1:] is ["cablam_file", "dssp_file", "Probe_file"  ## [1:] is a string slicing notation adding everything that follows the command itself.

args = sys.argv[1:]     

#sys.argv is a [] list of everything written on the commandline
#sys.argv[0] is the command itself (in this case "python Angle_add_HACO_02292016.py)
#  so sys.argv[1:] is a list of everything after the initial command;
#  PDB files can be run individually this way by adding the PDB file name as sysargv[1:]

if len(args) == 0:       #no file submitted, print help text
  print >> sys.stderr, """
This script accepts a probe-generated text file as input.
That file must start with the header line which contains  "name:pat:type:srcAtom: ..."" 
"""
  sys.exit()           # Test worked; printed the above when NO file was entered.

  #print >> printing location, thing to print
  #  an alternate printing syntax useful for going places other than sys.stdout
  #  The trailing comma after 'sys.stderr,' is necessary.

cablam_3CAangles_filename = args[0]             ## see line 7 above: the cablam file will be the FIRST [0] argument following the command 
cablamFile = open(cablam_3CAangles_filename)  ## cablamFile will open the file named in the FIRST argument shown on the command line    

CablamData = ["ResID", "MuIn", "MuOut", "VirtDihed", "CA3angle"]  ## These are the 5 comma-delimited segments of cablam text output

cablamDict = {}

for line in cablamFile:
  if line.startswith(",pdb:model:chain"):        
    continue
  adjustlength = line.split(",")      ## cablam output is strictly comma-delimited and can't be set up in a space-oriented format at this stage.
  ResID = adjustlength[0]
  CA3angle = adjustlength[4].strip()
  key = ResID[4:10]   ## had to put precise spaces to match dssp format exactly
  
      ##cablamDict[key] = round(float(CA3angle), 3)
  if CA3angle == "NULL":
  	continue
  
  cablamDict[key] = "%.1f" % float(CA3angle)   ## Have to use the dictionary locator, not the variable name

## print ResID + "  " + CA3angle

## cablamFile.close()    ## Tests ok to here. 

##keys = DsspDict.keys()
##keys.sort()
##for key in keys:
##  print key, DsspDict[key]    ## This tested perfectly: each ResNumber + ChainID is a key, and each SecStruct is a DsspDict key. 

##import sys
##sys.exit()


dssp_file_name = args[1]
dssp_file = open(dssp_file_name)

DsspData = ['H','G','I','E','B','S','T']

DsspDict = {} 

#Results = {'H':[],'G':[],'I':[],'E':[],'B':[],'S':[],'T':[]}
#this dictionary will hold the ResNumIDs pulled from the dssp file
#  it is keyed by the dssp letter codes

for line in dssp_file:
  ChainID = line[11:12]
  ResNumber = line[6:10]
  ResName = line[13:15]
  SecStruct = line[16:17]
  key = ChainID + ResNumber + " "
  
    
  DsspDict[key] = SecStruct     ##adds this to dict: {key:SecStruct}

  ##print key + "  " + cablamDict[key] + " " + line[16:17] 



probe_file_name = args[2]
probefile = open(probe_file_name)

ProbeData = ["PDBID", "contact", "SrcChainID", "SrcResNum", "SrcResName", "SrcAtom", "TgtChainID", "TgtResNum", "TgtResName", "TgtAtom", "GapWidth" ]

ProbeDict = {}

#Results = {'H':[],'G':[],'I':[],'E':[],'B':[],'S':[],'T':[]}                                                                                            
#this dictionary will hold the ResNumIDs pulled from the dssp file                                                                                       
#  it is keyed by the dssp letter codes                                                                                                                  

for line in probefile:
  if line.startswith("name"):
    continue
  adjustlength = line.split(":")
  SrcResNum = line[11:15]     
  SrcChainID = line[10:11]
  SrcResName = line[16:19]
  SrcAtom = line[20:25]
  TgtResNum = line[28:32]
  TgtChainID = line[27:28]
  TgtResName = line[33:36]
  TgtAtom = line[37:41]
  GapWidth = adjustlength[6]
  Srckey = SrcChainID + SrcResNum + " "
  Tgtkey = TgtChainID + TgtResNum + " "

  ProbeDict[Srckey] = GapWidth          ##adds this to dict: {key:SecStruct}                                                                               \
                                                                                                                                                         
  ## print ProbeDict[key]
  
  if Srckey not in cablamDict:
    continue
  if Tgtkey not in cablamDict:
    continue
    
  ##keys = DsspDict.keys()
##keys.sort()
##for key in keys:
##  print key, DsspDict[key]  


  print "1hm9H" + "," + Srckey + "," + SrcResName + "," + SrcAtom + "," + cablamDict[Srckey] + "," + DsspDict[Srckey] + "," + TgtChainID + "," + TgtResNum + "," +  TgtResName + "," + TgtAtom + "," + cablamDict[Tgtkey] + "," + DsspDict[Tgtkey] + "," + ProbeDict[Srckey]

probefile.close()
dssp_file.close()
cablamFile.close()

##  START HERE - figure out how to add back in the DSSP designations and add Cablam 3CA angle for Tgt residues .  Mark on excel sheets all crease residues and all beta-barrel dots residues.
