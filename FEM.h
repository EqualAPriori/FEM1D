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
	}
	void printVec(vector<double>& vec){
		int ii;
		for(ii=0; ii<vec.size(); ii++){
			cout << vec[ii] << " ";
		}
		cout << endl;
	}
	void printVec2D(vector<vector<double> >& mat){
		int irow, icol;
		for(irow=0; irow<mat.size(); irow++){
			for(icol=0; icol<mat[irow].size(); icol++){
				cout << mat[irow][icol] << " ";
			}
			cout << endl;
		}
	}
	void printVec2D(vector<vector<int> >& mat){
		int irow, icol;
		for(irow=0; irow<mat.size(); irow++){
			for(icol=0; icol<mat[irow].size(); icol++){
				cout << mat[irow][icol] << " ";
			}
			cout << endl;
		}
	}

//public:
	FEM(int ord);
	//~FEM();

	// Initializations
	void setGrid(vector<double>& zcoords);
	void GetGLL(int npts, vector<double>* xi, vector<double>* w, vector<vector<double> >* h); //write results into {xi,w,h}
	
	// Building and Assembling calculation
	void BuildLaplacian(vector<double>& eps, vector<double>* D2out);
	void BuildLaplacian(vector<double>& eps){ BuildLaplacian(eps, &D2); };
	void BuildLaplacian(vector<double>* D2out);
	void BuildLaplacian(){ BuildLaplacian(&D2); };
	
	void BuildSource(vector<double>& source, vector<double>* fout);
	void BuildSource(vector<double>& source){ BuildSource(source, &f); }
	void BuildMass(vector<double>& mass, vector<double>* Mout);
	void BuildMass(vector<double>& mass){ BuildMass(mass, &M); }

	void BuildBoundary(double const bc[2][3], vector<double>* A, vector<double>* b);
	void BuildBoundary(double const bc[2][3]){ BuildBoundary(bc, &Op, &f); };
	void BuildBoundary(){ BuildBoundary(BCs, &Op, &f); }

	void BuildIntegral(vector<double>& kernel, vector<double> *Kout);
	void BuildIntegral(vector<double>& kernel){ BuildIntegral(kernel,&K); }

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



