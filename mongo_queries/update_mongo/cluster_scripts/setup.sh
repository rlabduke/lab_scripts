#! /bin/bash

# make the requisite directories for the output
mkdir val_out
mkdir results

# write the correct full path is qsub.sh
sed -i 's:REPLACEME:'`pwd`':' qsub.sh

# get the pdbs that need to be done 
python get_missing.py residues_colkeys_pdbs.l

# write out the number of pdbs to be done
echo
echo "This is how many lines are in pdbs_to_be_done.l"
wc -l pdbs_to_be_done.l
echo
echo

a=($(wc pdbs_to_be_done.l))
lines=${a[0]}
words=${a[1]}
chars=${a[2]}
echo "Next do this (make sure phenix is sourced first):"
echo "qsub -t 1-$lines qsub.sh"
