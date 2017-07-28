#include <iostream>
#include "FEM.h"
//#include <Accelerate/Accelerate.h>
using namespace std;

extern "C" void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);


int main(){
    FEM f(2);
    cout << "size " << f.ws.size() << " test " << f.xis.at(1) << " weight " << f.ws.at(1) << endl;
    cout << "h " << f.hs[0][0] << ' ' << f.hs[0][1] << ' ' << f.hs[1][0] << ' ' << f.hs[1][1] << endl;
    
    cout << endl << "Grid Initialization Test" << endl;

    vector<double> zs = {0,1,5};//{0,1,2,3,4,5};
    f.setGrid(zs);
	f.printVec(f.zs);

	cout << endl;

	f.printVec2D(f.indOfElemE);

    cout << endl;

    f.printVec2D(f.zsOfElemE);

    cout << endl << "Laplacian Construction Test" << endl;
    //double eps[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    //zs.assign(eps,eps+16);
    vector<double> eps(f.Nz,1);
    zs.assign(eps.begin(),eps.end());
    f.BuildLaplacian(zs);

    f.printSqMat(f.D2);

    cout << endl<< "Source Test" << endl;
    
    fill(eps.begin(),eps.end(),1);
    f.BuildSource(eps,&f.f);
    f.printVec(f.f);

    cout << endl<< "Mass Test" << endl;

    f.printVec(eps);
    f.BuildMass(eps,&f.M);
    f.printVec(f.M);

    cout << endl << "Boundary Initialization Test" << endl;

    f.Op.assign(f.D2.begin(),f.D2.end());
    f.BuildBoundary();
    f.printSqMat(f.Op);
    cout<<endl;
    f.printVec(f.f);

    cout << endl << "Integral Kernel Construction Test" << endl;
    vector<double> k(f.Nz*f.Nz,1);
    f.BuildIntegral(k);
    f.printSqMat(f.K);


    return 0;
}

