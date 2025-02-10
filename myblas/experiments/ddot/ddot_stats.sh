#!/bin/bash

if (($#==0))
    then
        Output=ddot_logs
else
    Output=$1    
fi
mkdir $Output
cd $Output

for((N=100;N<=500;N+=100))
do
    echo "ddot with N=" $N 
    sbatch --job-name=ddot_N${N} \
    --error=../err_${Output}.log --open-mode=append ../ddot_job.slm $N $(cd ../../../../build/default/testings/ && pwd)
done

