# Custom BLAS


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