#!/bin/bash

#/home/smlewis/whole_PDB_as_PDB_20160429_gunzip is directory with whole PDB, gunzipped

pdbpath=/home/smlewis/whole_PDB_as_PDB_20160429_gunzip

find $pdbpath -name "*ent" -exec grep -H ^SSBOND {} \; > all_ssbond_recs.txt

cat all_ssbond_recs.txt | sed 's>'"$pdbpath"'/[0-9a-z][0-9a-z]/pdb>>g' | sed 's>.ent:>.pdb >g' > all_ssbond_recs.proc.txt

