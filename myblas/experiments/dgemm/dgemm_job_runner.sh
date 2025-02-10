#!/bin/bash

if (($#==0))
    then
        echo "Argument manquant (perf ou warm)"
        exit
else
    cache=$1
fi

#[ ! -d "dgemm_log"] && mkdir dgemm_log

mkdir dgemm_log

cd dgemm_log

sbatch --job-name=dgemm --error=../err_dgemm.log --open-mode=append ../dgemm_exp_job.slm $cache


