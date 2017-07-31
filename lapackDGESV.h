// Class object for interfacing use of lapack's DGESV matrix solves.
// Right now defaults for handling only NRHS = 1;
// Also usign default values for various flags and options.

#ifndef __DGESV__
#define __DGESV__

#include <cmath>
#include <iostream>
#include <vector>

extern "C" void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);

extern "C" void dgesvx_(char* FACT, char* TRANS, int* N, int* NHRS, double* A, int* LDA, double* AF, int* LDAF,\
		int* IPIV, char* EQED, double* R, double* C, double* B, int* LDB, double* X, int* LDX,\
		double* RCOND, double* FERR, double* BERR, double* WORK, int* IWORK, int* INFO);

using namespace std;

class lapackDGESV{
public:
	int N, NRHS, LDA, LDB, INFO, LDAF, LDX;
	int *IPIV, *IWORK;
	double RCOND;
	double *AF, *R, *C, *FERR, *BERR, *WORK;

	double *A, *B, *X; 	//these will be fed in as arguments to the solver function.
						//not sure if I want this object to keep a record of where its answers are...
	char FACT = 'N';
	char TRANS = 'N';
	char EQUED = 'N';

	bool allocated = 0;

	//METHODS:
	lapackDGESV():NRHS(1){}
	lapackDGESV(int nn, int nnrhs=1){
		refresh(nn, nnrhs);
	}
	void refresh(int nn, int nnrhs=1){
		N = nn;
		NRHS = nnrhs;
		LDA = N;
		LDB = N;
		LDAF = N;
		LDX = N;

		if(allocated == 1){clear();}
		IPIV = new int[N];
		IWORK = new int[N];

		AF = new double[N*N];
		R = new double[N];
		C = new double[N];
		FERR = new double[NRHS];
		BERR = new double[NRHS];
		WORK = new double[4*N];

		allocated = 1;
	}
	void clear(){
		delete [] IPIV;
		delete [] IWORK;

		delete [] AF;
		delete [] R;
		delete [] C;

		delete [] FERR;
		delete [] BERR;
		delete [] WORK;
	}

	//solve in place using dgesv_
	void solve(double *a, double *b){
		A = a;
		B = b;
		dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
	}
	//solve without overwriting
	void solve(double *a, double*b, double *x){
		A = a;
		B = b;
		X = x;
		dgesvx_(&FACT, &TRANS, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, &EQUED, R, C, B, &LDB,\
			X, &LDX, &RCOND, FERR, BERR, WORK, IWORK, &INFO );
	}

	~lapackDGESV(){
		clear();
	}
};
#endif






