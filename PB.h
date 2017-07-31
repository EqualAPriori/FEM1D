// Header for PB class
// default is 1:1 electrolyte, lb=10/4/pi nm, neumann on left bdry, dirichlet on right bdry
// 0.05/nm^3 conc, box is 0 to 10nm, sigma=0.01, no image charge.
//
//
#ifndef __PB__
#define __PB__

#include "FEM.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "lapackDGESV.h"

// This function call MODIFIES A: 
//	On exit, the factors L and U from the factorization A = P*L*U; 
//	the unit diagonal elements of L are not stored.
extern "C" void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
// This iterative solver had trouble solving... and hence doesn't leave matrix untouched
// extern "C" void dgesv_(int* N, int* NHRS, double* A, int* LDA, int* IPIV, double* B, int* LDB, double* X, int* LDX, double* WORK, float* SWORK, int* iter, int* INFO);
// Another moethod
extern "C" void dgesvx_(char* FACT, char* TRANS, int* N, int* NHRS, double* A, int* LDA, double* AF, int* LDAF,\
		int* IPIV, char* EQED, double* R, double* C, double* B, int* LDB, double* X, int* LDX,\
		double* RCOND, double* FERR, double* BERR, double* WORK, int* IWORK, int* INFO);


using namespace std;

class PB{
public:
//private:
	FEM fpb; 						//default is 2nd order in this PB class (see constructor)

	int Nspecies = 2;
	vector<double> q {1,-1}; 		// valence

	vector<vector<double> > rho;
	vector<double> rho0{0.05,0.05}; // default conc.
	vector<double> netchg;
	vector<double> netchgLessDH; 		//subtracts off the DH term and scales RHS by 4*pi*lb (since I don't have epsilon in Laplacian)
	vector<double> kap2;
	double kap20;
	double lb=10.0/4.0/M_PI;
	double sigma = 0.01;			// default surface charge.
	double BCs[2][3] = {{1,0,-0.1},{0,0,0}}; //default

	int order = 2;
	double L = 10.0;
	int Nz = 101;
	double dx = 0.1;

	vector<double> nodes;			// coarse node list. FEM may refine it.
	vector<double> psi;

	int iterMax = 1e3; 
	int iterations = 0;
	double tol = 1e-10; 
	double mixBase = 0.1; 			// "error "tolerance" for using mixing rate 1"
	double mix = 1.0;
	double err = 1.0;
	int errflag = 0;
	lapackDGESV myDGESV;


//public:
	PB();
	PB(int ord);
	//~PB();

	void setGrid(){ setGrid(L,Nz); }
	void setGrid(double Len, int Nx);
	void setGrid(vector<double>& nodes);
	void setBC(double const bc[2][3]);

	void bulkInit();
	void updateDensity();
	void setOperator();
	void solve1step();
	void solve();

	void nextSol(vector<double>& psiNew, int flag=0); //either Picard or DIIS or something else
	double calcError(vector<double>& ynew, vector<double>& yold);
};

#endif //defined __PB__


