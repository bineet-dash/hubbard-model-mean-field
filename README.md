# MC Sampling for Half-filling Hubbard Model.
This is a C++ implementation of Monte-Carlo sampling for SPA Hamiltonian of a Half-filling Hubbard Model. The project uses Eigen Library for Matrix Manipulation.

# System Requirements:
 1. C++ (g++ -std=c++14)
 2. Eigen Library
 3. Lapack Library (llapack, llapacke)

# Installation Instruction for Eigen
Please download Eigen tarball from http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2 . Follow the [instruction manual](http://eigen.tuxfamily.org/dox/GettingStarted.html#title2) regarding installation of Eigen.
(**Linux and Mac OS X users:** download and extract the tar ball, then create a symlink in /usr/local/include to the "Eigen" folder.)

Install Lapack library, if not already installed. Most Linux distributions have made Lapack available their official repositories. You can follow [this guide](http://icl.cs.utk.edu/lapack-for-windows/lapack/) for installation in Windows. 
(**Ubuntu users** can install with apt-get:  _sudo apt-get install liblapack-dev liblapack3 liblapacke liblapacke-dev_) 

To build the program compile using "g++ -std=c++14 _corresponding source file_ -llapack -llapacke".
