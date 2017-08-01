// Header for Green class
// Note that Green ftn doesn't inherently know anything about the electrolyte solution,
// which enters only through the ionic strength function
// Green's ftn only needs to know about the dielectric constant and the B.C.

// default calculation right now is one wall, with dirichlet at right B.C. 
// (maybe better to use slope going to zero... implicitly becomes two-wall problem?)
// see Rui Wang JCP 2015 for B.C.


#ifndef __Green__
#define __Green__

#include "FEM.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "lapackDGESV.h"

using namespace std;

class Green{
public:
//private:
	FEM fe; 						//default is 2nd order in this Greenftn class (see constructor)

	double lb=10.0/4.0/M_PI;
	double lbwall = 0.0;
	double BCs[2][3] = {{1,0,0},{0,0,0}}; //default

	int order = 2;
	double L = 20.0;
	int Nz = 101;
	double dx = 0.1;

	double k2prev = 0.0;			// keeping track for easy undoing later on. Note is k^2
	vector<double> diagKernelPrev;
	vector<double> nodes;			// coarse node list. FEM may refine it.
	vector<double> psi;				// the most recently solved instance of the Green's ftn

	int errflag = 0;
	lapackDGESV myDGESV;

//public:
	Green();
	Green(int ord);
	//~Green();

	void setGrid(){ setGrid(L,Nz); }
	void setGrid(double Len, int Nx);
	void setGrid(vector<double>& nodes);
	void setBC(double const bc[2][3]);
	void setlb(double lbs, double lbw);

	// note that this `eps` is probably going to be a scaling relative to the solution lb,
	// which is how my PB function is implemented, where I pre-multiplied the 4*pi*lb to the RHS 
	// to turn density --> kappa^2
	void BuildBareOperator(vector<double> &kernel, double eps=1.0);
	void BuildBareOperator(vector<double> &kernel, vector<double> &eps);
	void modifyDeltaFtnKernel(vector<double> &kernel, vector<double> *Op);
	void modifyDeltaFtnKernel(vector<double> &kernel){ modifyDeltaFtnKernel(kernel, &fe.Op); }
	void modifyDeltaFtnKernel(double kernel, vector<double> *Op);
	void modifyDeltaFtnKernel(double kernel){ modifyDeltaFtnKernel(kernel, &fe.Op); }
	void kOperator(double k, vector<double> *Op);
	void kOperator(double k){ kOperator(k, &fe.Op); }

	void solve(double z0, double k);
	void solve(int index, double k);

	void buildDelta(int loc, vector<double> *f);
	void buildDelta(int loc){ buildDelta(loc, &fe.f);}
	void buildDelta(double z0, vector<double> *f);
	void buildDelta(double z0){ buildDelta(z0, &fe.f); }

	double find(double z0, vector<double> &pos);
	double find(double z0){ return find(z0, fe.zs); }
};

#endif //defined __PB__


