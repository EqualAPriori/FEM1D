# FEM1D


Compile with (Mac):
clang++ -std=c++11 -L/usr/local/opt/openblas/lib -lopenblas -llapack  -O3 main.cpp FEM.cpp -o ./test


Note currently coded s.t. a FEM class instance can be built to build matrices/write results to arrays whose addresses are pointed in (opearating only on object's own data). However, this isn't robust/complete yet because the methods still implicitly assume that the matrices fed in are commensurate with the spatial grid stored in the FEM object. "Normal" class behavior of modifying object's own data is provided via overloaded constructors although if feeding in user matrices, but one needs to be careful that they are 1) commensurate in size with the spatial grid in the object instane, and 2) properly zeroed first since, where possible, I only implement sparse zeroing.



