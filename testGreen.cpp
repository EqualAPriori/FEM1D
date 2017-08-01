#include <iostream>
#include <iomanip>
#include <fstream>
#include "lapackDGESV.h"
#include "PB.h"
#include "Green.h"
//#include <Accelerate/Accelerate.h>
using namespace std;



int main(){
	/* Testing my BuildIntDiffOperator function
	PB mypb(1);
	vector<double> z{0,1,2,3};
	mypb.setGrid(z);
	vector<double> k{1,0,0,0, 0,2,0,0, 0,0,3,0, 0,0,0,4};
	mypb.fpb.BuildIntegral(k);
	mypb.fpb.printSqMat(mypb.fpb.K);
	mypb.fpb.BuildIntDiffOperator(k,1,0.1);
	mypb.fpb.printSqMat(mypb.fpb.Op);
	//*/

	// Checking the internals
	vector<double> z{0,1,2,3,4};
	Green g(1);
	g.setGrid(z);
	g.buildDelta(2.1);
	g.fe.printVec(g.fe.zs);
	g.fe.printVec(g.fe.f);

	vector<double> k(z.size(),0.1);
	//vector<double> k = {1,0,0,0,0, 0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1};
	g.BuildBareOperator(k);
	g.fe.printSqMat(g.fe.Op);
	g.kOperator(1);
	g.fe.printSqMat(g.fe.Op);
	g.kOperator(1);
	g.fe.printSqMat(g.fe.Op);
	g.kOperator(10);
	g.fe.printSqMat(g.fe.Op);


	// Try solving a green's function calculation
	vector<double> z1;
	z1.reserve(601);
	double x;
	for(x=0.0; x<50.0; x+=0.5){
		z1.push_back(x);
	}
	Green g1(2);
	g1.setGrid(z1);
	vector<double> kap2(g1.fe.Nz,0.1);
	g1.BuildBareOperator(kap2);
	double z0 = 25.0;
	double kw = sqrt(0.1);
	g1.solve(z0, kw);
	g1.fe.printVec(g1.psi);


	///*
	ofstream outFile;
    outFile<<setiosflags(ios::left)<<setiosflags(ios::fixed);
    outFile.open("out.txt");

    int ii;
    for(ii=0; ii<g1.fe.Nz; ii++){
    	outFile << g1.fe.zs[ii] << "\t" << g1.psi[ii] << endl;
    }

    outFile.close();
    outFile<<resetiosflags(ios::left)<<resetiosflags(ios::fixed);
    //*/


    return 0;
}