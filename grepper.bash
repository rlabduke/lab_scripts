#!/bin/bash

#/home/smlewis/whole_PDB_as_PDB_20160429_gunzip is directory with whole PDB, gunzipped

pdbpath=/home/smlewis/whole_PDB_as_PDB_20160429_gunzip

find $pdbpath -name "*ent" -exec grep -H ^SSBOND {} \; > all_ssbond_recs.txt
