# MC Sampling for Half-filling Hubbard Model.
This is a C++ implementation of Monte-Carlo sampling for SPA Hamiltonian of a Half-filling Hubbard Model. The project uses Eigen Library for Matrix Manipulation.


# System Requirements:
 1. C++ (g++ -std=c++14)
 2. Eigen Library


# Installation Instruction for Eigen
Please download Eigen tarball from http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2 . Follow the instruction manual from http://eigen.tuxfamily.org/dox/GettingStarted.html#title2 regarding installation of Eigen.


(**TLDR for Linux and Mac OS X users:** download and extract the tar ball, then create a symlink in /usr/local/include to the "Eigen" folder.)

To build the program compile using "g++ -std=c++14 main.cpp spalib.cpp".
To enable debugging functions,include "debug.h" and compile with "g++ -std=c++14 main.cpp spalib.cpp debug.cpp".
