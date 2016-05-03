#!/usr/bin/python

'''This script parses SSBOND records and attempts to find N to N+1 disulfide bonds.'''

#expected record format:
#1isz.pdb SSBOND   2 CYS A  254    CYS A  260                          1555   1555  2.03  
#Note this is fixed width
#Note that an insertion code MAY occur, so we are processing as fixed-width

#read file (from arg?  hardcode)
filename = "all_ssbond_recs.proc.txt"
with open(filename) as f:
    for SSBOND in f:
        #just skip if insertion codes involved
        space = " "
        icode1 = SSBOND[30]
        icode2 = SSBOND[44]
        if (icode1 != space) or (icode2 != space):
            #print icode1, "icodespacer", icode2
            continue
        chain1=SSBOND[24]
        chain2=SSBOND[38]
        #print chain1, "chainspacer", chain2
        if chain1 == chain2:
            #these counts would include insertion code
            #resn1=SSBOND[26:31]
            #resn2=SSBOND[40:45]
            resn1=int(SSBOND[26:30])
            resn2=int(SSBOND[40:44])
            #print resn1, "resnspacer", resn2
            if (resn1+1) == resn2:
                print SSBOND.rstrip()
                #pass


#remove hits from other list
#delete anything with an insertion code
