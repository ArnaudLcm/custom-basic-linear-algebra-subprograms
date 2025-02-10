#!/bin/bash

wd=$(cd ../../../../build/default/testings/ && pwd)
cache=$1

guix shell --pure -D mini-chameleon --with-input=openblas=mkl bash gcc-toolchain emacs nano vim --tune -- bash --norc
cd $wd


# faire nos propres benchmarks pour tester la variation de taille de blocs ?

N=1000
BMAX=1200



for((n=100;n<$N;n+=100))
do
for b in {1,2,4,5,10}
do


            BLOCKSIZE=$(($n/$b)) ./${cache}_dgemm -v seq -N $n
            
            echo -e
    

done


done


N=2000
BMAX=1200



for((n=1000;n<$N;n+=200))
do
for b in {1,2,4,5,10}
do


            BLOCKSIZE=$(($n/$b)) ./${cache}_dgemm -v seq -N $n
            
            echo -e
    

done


done

N=3000
BMAX=1200



for((n=2000;n<$N;n+=300))
do
for b in {1,2,4,5,10}
do


            BLOCKSIZE=$(($n/$b)) ./${cache}_dgemm -v seq -N $n
            
            echo -e
    

done


done
