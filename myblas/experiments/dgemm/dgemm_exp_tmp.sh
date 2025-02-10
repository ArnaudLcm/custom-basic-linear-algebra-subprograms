#!/bin/bash

wd=$(cd ../../../../build/default/testings/ && pwd)
cache=$1

guix shell --pure -D mini-chameleon --with-input=openblas=mkl bash gcc-toolchain emacs nano vim --tune -- bash --norc
cd $wd


# faire nos propres benchmarks pour tester la variation de taille de blocs ?

N=3000
S=0

for b in {3600,3000}
do
for((k=0;k<6;k++))
do

  
            ./${cache}_dgemm -v seq -N $b
            echo -e

        
done
done