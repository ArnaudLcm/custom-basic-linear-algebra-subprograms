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

for((n=50;n<1200;n+=n/4))
do
    BLOCKSIZE=$n ./${cache}_dgemm -v seq -N $n
    echo -e
    S=$n
done


for((n=$S+$S/4;n<$N;n+=n/4))
do
for((i=1250;i>0;i--))
do
         if(($n%$i==0))
         then   
            BLOCKSIZE=$i ./${cache}_dgemm -v seq -N $n
            echo -e
            break
        fi
done
done
done