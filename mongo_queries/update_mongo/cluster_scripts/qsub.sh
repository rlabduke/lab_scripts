#! /bin/csh -q
#$ -cwd
#$ -o residueval3.output -j y -N residueval
limit datasize 2000000

source /net/chevy/raid1/bhintze/phenix_clean/build/setpaths.csh



phenix.python REPLACEME/qsub.py $SGE_TASK_ID $SGE_TASK_LAST >& results/residueval.$SGE_TASK_ID.out



exit

