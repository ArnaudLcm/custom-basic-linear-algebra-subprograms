#!/bin/bash

if (($#==0))
    then
        echo "Argument manquant (perf ou warm)"
        exit
else
    cache=$1
fi

mkdir dgetrf_log
cd dgetrf_log

sbatch --job-name=dgetrf --error=../err_dgetrf.log --open-mode=append ../dgetrf_exp_job.slm $cache


