#!/bin/bash

if (($#==0))
    then
        echo "Argument manquant (perf ou warm)"
        exit
else
    cache=$1
fi

mkdir ddot_log
cd ddot_log

sbatch --job-name=ddot --error=../err_ddot.log --open-mode=append ../ddot_exp_job.slm $cache


