# Mantle_Convection_Lab
This is a 2D mantle convection code, written in C with finite element method, parallel computing, and adaptive mesh. This code is written by Mingming Li at Arizona State University.
This is just the first version. There remains a lot to be added to this code. The goal for this code is to deal with multi-scale dynamics and complex rheology more effectively. The code is self-consistent and does not rely on any other software.


To compile the code, install mpi, and type 'make' on the command line.
To run the test model, type 'mpirun -np 16 ./MCC.mpi mcc01.input'. The results will be saved under the folder of mcc01
