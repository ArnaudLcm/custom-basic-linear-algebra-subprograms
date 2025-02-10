#!/bin/bash

wd=$(cd ../../../../build/default/testings/ && pwd)
cache=$1





guix shell --pure -D mini-chameleon --with-input=openblas=mkl bash gcc-toolchain emacs nano vim --tune -- bash --norc
cd $wd



# N=1000


# for((b=50;b<$N;b+=50))
# do
#     r=$N%$b
#     if (($r==0))
#         then
#             ./${cache}_dgemm -v seq -N 1000 -b --nb=$b
#             echo -e
#     fi
# done

# for((N=50;N<3000;N+=50))
# do 
#     r=$N%4
#     if (($r==0))
#         then 
#             ./${cache}_dgemm -v seq -N $N
#             echo -e
#     fi
# done

# echo $#

nb_threads=$2
for ((N=50;N<3000;N+=50))
do 
    OMP_NUM_THREADS=$nb_threads ./${cache}_dgemm -v omp -N $N
    echo -e
done

# if (($#==2))
#     then
#         for ((N=50;N<3000;N+=50))
#         do 
#             ./${cache}_dgemm -v seq -N $N
#             echo -e
#         done
# fi





