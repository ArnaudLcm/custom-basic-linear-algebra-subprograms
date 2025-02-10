#!/bin/bash

wd=$(cd ../../../../build/default/testings/ && pwd)
cache=$1

guix shell --pure -D mini-chameleon --with-input=openblas=mkl bash gcc-toolchain emacs nano vim --tune -- bash --norc
cd $wd


# faire nos propres benchmarks pour tester la variation de taille de blocs ?

N=5000
S=0

for((h = 0;h<6;h++))
do

for((n=64;n<1200;n=n*2))
do
    BLOCKSIZE=$n ./${cache}_dgemm -v seq -N $n
    echo -e
    S=$n
done


for((n=$S*2;n<$N;n=n*2))
do


            BLOCKSIZE=1024 ./${cache}_dgemm -v seq -N $n
            echo -e

done
done