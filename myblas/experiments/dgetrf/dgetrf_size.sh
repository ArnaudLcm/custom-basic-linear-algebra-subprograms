#!/bin/bash

if (($#==0))
    then
        Output=dgetrf_log
else
    Output=$1    
fi
mkdir $Output
cd $Output

# for((N=50;N<=5000;N+=N/4))
# do
#     for((M=50;M<=50000;M+=M/4))
#     do
#         for((K=50;K<=5000;K+=K/4))
#         do
#             echo "dgemm M=" $M " N=" $N " K=" $K 
#             sbatch --job-name=dgemm_job \
#             --error=../err_${Output}.log --open-mode=append ../dgemm_job.slm $M $N $K $(cd ../../../build/default/testings/ && pwd)
#         done
#     done
# done

for((N=50;N<=1000;N+=50))
do
    echo "dgemm M=" $N " N=" $N " K=" $N
    sbatch --job-name=dgetrf_job \
    --error=../err_${Output}.log --open-mode=append ../dgetrf_job.slm $N $N $N $(cd ../../../../build/default/testings/ && pwd)
done

