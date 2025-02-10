# Custom BLAS
## Introduction
Project made at the university to build a custom basic linear algebra library.

**Kernels available**:
- Sequential
- Distributed around cores through OpenMP
- Distributed around cpus through MPI
- Accelerated through a GPU

**Operations implemented**:
- DDOT
- DGEMM: Matrix operation C  :=alpha*op( A )*op( B ) + beta*C
- DGETRF: Computes an LU factorization of a general M-by-N matrix A with A = P * L * U



## Setup
Please verify that you have the following packages installed :
- CMake
- Guix

In order to install other depedencies (such as OpenMPI), Guix will be used.

First, please add this symlink to your guix config folder :
```bash
ln -s ./channels.scm ~/.config/guix #If you follow the XDG standard, otherwise your guix config folder
```
You can now pull depedencies from guix :
```bash
guix pull --allow-downgrades
hash guix #Make sure guix has the correct version
```


## Build

First, please create this folder : 
```bash
mkdir -p build/default
```

Launch the guix environment :
```bash
guix shell --pure -D mini-chameleon --with-input=openblas=mkl bash emacs nano vim -- bash --norc
```


Then, you can build the project :
```bash
cmake . -B build/default -DENABLE_MPI=ON -DENABLE_STARPU=ON
cmake --build build/default
```