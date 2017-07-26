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

	//GLL info
	vector<double> xis;
	vector<double> ws;
	vector<vector<double> > hs;

	//Calculation matrices and vectors -- should implement with simple arrays instead...
	vector<vector<double> > D2loc;	//temp variable, local to element
	vector<vector<double> > D2; 		//Laplacian
	vector<vector<double> > M;  		//Mass
	vector<vector<double> > K;  		//Integration Kernel
	vector<vector<double> > K2; 
	vector<vector<double> > II; 


//public:
	FEM(int ord);
	//~FEM();

	void setGrid(vector<double>& zcoords);
	void GetGLL(int npts, vector<double>* xi, vector<double>* w, vector<vector<double> >* h); //write results into {xi,w,h}
	
	//interestingly, if returning vector<double> can't separate into header and implementation...
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



//Ways to use data with variable size: https://stackoverflow.com/questions/751878/determine-array-size-in-constructor-initializer
//http://www.cplusplus.com/forum/articles/7459/
//but even vectors benefit from initializing size: https://stackoverflow.com/questions/25937413/assigning-values-to-2d-vector-by-using-indices
//https://www.quora.com/How-do-I-create-a-2D-vector-array-using-std-vector
/*

--- Templating your own class:

template < int ARRAY_LEN > // you can even set to a default value here of C++'11

class MyClass
  {   
  int array[ARRAY_LEN]; // Don't need to alloc or dealloc in structure!  Works like you imagine!   
  }
// Then you set the length of each object where you declare the object, e.g.
MyClass<1024> instance; // But only works for constant values, i.e. known to compiler

--- Being careful to allocate memory:

class MyClass
  {
  int *array;

  MyClass(int len) { array = calloc(sizeof(int), len); assert(array); }   
  ~MyClass() { free(array); array = NULL; } // DON'T FORGET TO FREE UP SPACE!
  }

--- Another option?

class Class
{
   int* array;
   Class(int x) : array(new int[x]) {};
};

*/