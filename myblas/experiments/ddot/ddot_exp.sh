#!/bin/bash

wd=$(cd ../../../../build/default/testings/ && pwd)
cache=$1

guix shell --pure -D mini-chameleon --with-input=openblas=mkl bash gcc-toolchain emacs nano vim --tune -- bash --norc
cd $wd

for((N=50;N<=500;N+=50))
do
    #echo "ddot with N=" $N 
    ./${cache}_ddot -v seq -N $N
    echo -e
done

