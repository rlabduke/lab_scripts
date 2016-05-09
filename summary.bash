#!/bin/bash

cat all_ssbond_recs.N-N+1pairs.txt | awk '{print $1}' | uniq > all_ssbond_recs.pdbs.txt
wc -l all_ssbond_recs.pdbs.txt
