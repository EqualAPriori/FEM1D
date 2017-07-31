// Header for FEM class FEM
#ifndef __FEM__
#define __FEM__
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

class FEM{
public:
//private:
	int order=1; 
	int NGLL = 2; //number of GLL points per element is 1+(order)

	int Nel, Nz;

	vector<double> zcoarse;
	vector<double> zs;
	vector<vector<double> > zsOfElemE;
	vector<vector<int> > indOfElemE;

	// GLL info
	vector<double> xis;
	vector<double> ws;
	vector<vector<double> > hs;

	// Calculation matrices and vectors -- should implement with simple arrays instead...
	// currently implemented in 1D column major form
	// to save memory, should contemplate building D2, M, K all directly onto single matrix variable
	double BCs[2][3] = {{0,1,0},{0,1,0}}; //default homogeneous Dirichlet
	vector<double> f;
	vector<vector<double> > D2loc;
	vector<double> D2;	// NxN
	vector<double> M;	// Nx1 since by construction diagonal with GLL quadrature evaluation
	vector<double> K;	// NxN
	vector<double> K2;	// NxN
	vector<double> II;	// NxN
	vector<double> Op; //the final assembled operator

	//vector<vector<double> > D2loc;	//temp variable, local to element
	//vector<vector<double> > D2; 		//Laplacian
	//vector<vector<double> > M;  		//Mass
	//vector<vector<double> > K;  		//Integration Kernel
	//vector<vector<double> > K2; 
	//vector<vector<double> > II; 


	// Column major indexing conversion and utilities
	int ij(int ii, int jj, int ncol){
		return ii + ncol*jj;
	}; 
	void printSqMat(vector<double>& mat){ //assuming is square
		int dim = sqrt(mat.size());
		int ii,jj;

		for(ii=0; ii<dim; ii++){
			for(jj=0; jj<dim; jj++){
				cout << mat[ij(ii,jj,dim)] << "\t";
			}
			cout << endl;
		}
		cout << endl;
	}
	void printVec(vector<double>& vec){
		int ii;
		for(ii=0; ii<vec.size(); ii++){
			cout << vec[ii] << " ";
		}
		cout << endl << endl;
	}
	void printVec2D(vector<vector<double> >& mat){
		int irow, icol;
		for(irow=0; irow<mat.size(); irow++){
			for(icol=0; icol<mat[irow].size(); icol++){
				cout << mat[irow][icol] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	void printVec2D(vector<vector<int> >& mat){
		int irow, icol;
		for(irow=0; irow<mat.size(); irow++){
			for(icol=0; icol<mat[irow].size(); icol++){
				cout << mat[irow][icol] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}


//public:
	FEM(int ord=1); 		//has default value s.t. can call empty
	void setOrder(int ord);
	//~FEM();

	// Initializations
	void setGrid(vector<double>& zcoords);
	void GetGLL(int npts, vector<double>* xi, vector<double>* w, vector<vector<double> >* h); //write results into {xi,w,h}
	
	// Building and Assembling calculation
	void BuildLaplacian(vector<double>& eps, vector<double>* D2out, bool clear=1);
	void BuildLaplacian(vector<double>& eps, bool clear=1){ BuildLaplacian(eps, &D2, clear); };
	void BuildLaplacian(vector<double>* D2out, bool clear=1, double eps=1.0 ); //unfortunately eps must be at end now.
	void BuildLaplacian(bool clear=1){ BuildLaplacian(&D2, clear); };
	
	void BuildSource(vector<double>& source, vector<double>* fout);
	void BuildSource(vector<double>& source){ BuildSource(source, &f); }
	void BuildSource(double source, vector<double>* fout);
	void BuildSource(double source){ BuildSource(source, &f); }

	void BuildMass(vector<double>& mass, vector<double>* Mout); //builds into Nx1 vector
	void BuildMass(vector<double>& mass){ BuildMass(mass, &M); }
	void BuildMass(double mass, vector<double> *Mout);
	void BuildMass(double mass){ BuildMass(mass, &M); }
	void BuildMassMatrix(vector<double>& mass, vector<double>* Mout); //builds into NxN matrix, *without* clearing
	void BuildMassMatrix(double mass, vector<double> *Mout); //builds into NxN matrix, without clearing

	void BuildBoundary(double const bc[2][3], vector<double>* A, vector<double>* b);
	void BuildBoundary(double const bc[2][3]){ BuildBoundary(bc, &Op, &f); };
	void BuildBoundary(){ BuildBoundary(BCs, &Op, &f); }
	void BuildBoundaryOperator(double const bc[2][3], vector<double>* A);
	void BuildBoundaryOperator(double const bc[2][3]){ BuildBoundaryOperator(bc, &Op); }
	void BuildBoundarySource(vector<double>* b);
	void BuildBoundarySource(){ BuildBoundarySource(&f); }

	void BuildIntegral(vector<double>& kernel, vector<double> *Kout);
	void BuildIntegral(vector<double>& kernel){ BuildIntegral(kernel,&K); }

	// Right now not gathering the integral kernel `K` yet.
	void combineOperator(vector<double>& Lap, vector<double>& mass, vector<double>* A);
	void combineOperator(){ combineOperator(D2, M, &Op); }

	// Updated version for Green's Function calculation to allow direct assembling
	// of operator into one matrix without zeroing or need for O(N^2) copying
	// Right now don't have version for scalar eps + vectorial mass.
	void BuildIntDiffOperator(vector<double> &kernel, vector<double> &eps, vector<double> &mass, vector<double> *OpOut);
	void BuildIntDiffOperator(vector<double> &kernel, vector<double> &eps, vector<double> &mass){
		BuildIntDiffOperator(kernel, eps, mass, &Op);
	};
	void BuildIntDiffOperator(vector<double> &kernel, vector<double> &eps, double mass, vector<double> *OpOut);
	void BuildIntDiffOperator(vector<double> &kernel, vector<double> &eps, double mass){
		BuildIntDiffOperator(kernel, eps, mass, &Op);
	}
	void BuildIntDiffOperator(vector<double> &kernel, double eps, vector<double> &mass, vector<double> *OpOut);
	void BuildIntDIffOperator(vector<double> &kernel, double eps, vector<double> &mass){
		BuildIntDiffOperator(kernel, eps, mass, &Op);
	}
	void BuildIntDiffOperator(vector<double> &kernel, double eps, double mass, vector<double> *OpOut);
	void BuildIntDiffOperator(vector<double> &kernel, double eps, double mass){
		BuildIntDiffOperator(kernel,eps,mass, &Op);
	}


	// Other utilities 
	// Interestingly, if returning vector<double> can't separate into header and implementation...
	vector<double> xiToZ(vector<double>& xi, double zleft, double zright){
		vector<double> z(xi.size());
		//z->assign(xi.begin(),xi.end());
		if(zleft >= zright){
			cout << "Error: zleft >= zright in xiToZ()" << endl;
			return z;
		}

		double len = fabs(zright - zleft);
		int ii;
		for(ii=0; ii < xi.size(); ii++){
			(z)[ii] = zleft + len*(xi[ii] + 1.0)/2.0;
		}
		return z;
	}
	void xiToZ(vector<double>& xi, double zleft, double zright, vector<double>* z);
	void zToXi(vector<double>& z, double zleft, double zright, vector<double>* xi);
};

#endif //defined (__FEM__)



