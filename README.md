# FEM1D
1D FEM code. Implemented a class that handles building the matrix operator, source, and boundary conditions for polynomial elements up to order 5. Higher order is also possible -- just need to code in the appropriate Gauss-Lobatto-Lgendre (GLL) nodes, weights, and derivatives (i.e. possibly read in from a file).

In addition, the FEM code has been used to solve new Poisson Boltzmann and electrostatic fluctuation (self energy) problems, based on the work of [Wang 2010, PRE 81 (021401)](http://link.aps.org/doi/10.1103/PhysRevE.81.021501) and [Wang and Wang 2015, JCP 142 (104705)](http://dx.doi.org/10.1063/1.4914170), where the Green's function (fluctuations) are the result of the ionic screening of point charges. We now do calculations within the same renormalied Gaussian fluctuation theory, but where the ionic screening can come from objects of finite size and shape. The importance of these non-point-charge screening objects was discussed for polyelectrolytes in [Shen and Wang 2017, JCP 146 (084901)](http://dx.doi.org/10.1063/1.4975777).

## Getting Started
Compilation command on (Mac):
clang++ -std=c++11 -L/usr/local/opt/openblas/lib -lopenblas -llapack  -O3 FEM.cpp PB.cpp Green.cpp main.cpp -o ./out

### Prerequisites:
This code requires the lapack and blas libraries, which are likely already included on many computers. You may use your own llapack or other blas libraries (we use openblas). In the future, more calculations may involve the Eigen or Armadillo libraries for ease of manipulating lots of matrix objects.

```
examples
```
and more
```
examples
```

The testPB.cpp, testGreen.cpp, and testFEM.cpp code give (hard-coded) examples of how to use the various objects, with terminal outputs to demonstrate what the matrices and vectors look like.

## Notes on Usage
* Note currently coded s.t. a FEM class instance can be built to build matrices/write results to arrays whose addresses are pointed in (opearating only on object's own data). However, this isn't robust/complete yet because the methods still implicitly assume that the matrices fed in are commensurate with the spatial grid stored in the FEM object. "Normal" class behavior of modifying object's own data is provided via overloaded constructors although if feeding in user matrices, but one needs to be careful that they are 1) commensurate in size with the spatial grid in the object instane, and 2) properly zeroed first since, where possible, I only implement sparse zeroing.


Objects:
* PB object -- Currently does primitive Poisson Boltzmann calculation. Can handle > 2 species
* Green's function object -- Calculates the Green's function given a scaled Ionic strength kernel, at wavenumber `k` and position `z0`
* FEM object -- handles building FEM matrices, vectors, and boundary conditions, given a polynomial `order` and coarse mesh/grid `zcoarse`

## Future
* Using all the classes in a full self-consistent self energy calculation
* A parser to read in parameters from a text file instead of hard-coding in.

## Contributing

## Versioning

## Authors
Kevin Shen

## Acknowledgments
* README template borrowed from [PurpleBooth](ttps://gist.github.com/PurpleBooth/109311bb0361f32d87a2#file-readme-template-md)
* FEM method referenced heavily from [Professor Jean Paul Ampuero](http://web.gps.caltech.edu/~ampuero/) and his [lecture notes](https://pdfs.semanticscholar.org/de95/5fb1f0d63f378f705ae8d97c9422f37e62f9.pdf)
