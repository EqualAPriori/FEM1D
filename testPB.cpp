#include <iostream>
#include <iomanip>
#include <fstream>
#include "lapackDGESV.h"
#include "PB.h"
//#include <Accelerate/Accelerate.h>
using namespace std;



int main(){
	PB mypb;

	//vector<double> z{0,0.00001,0.0005,0.0025,0.01,0.025,0.05,0.1,0.2,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,10.0};
	vector<double> z{0,0.01,0.025,0.05,0.1,0.2,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,10.0};
	//vector<double> z{0,0.5,1.0,1.5,2.0,2.5,3.0};
	//mypb.setGrid();
	//mypb.setGrid(2,11);
	double BCs[2][3] = {{1,0,-0.1},{0,0,0}}; //not very happy with stronger surface charge... optimal mixing ~0.7 for sfchg = -0.3
	mypb.setBC(BCs);

	mypb.setGrid(z);
	mypb.bulkInit();
	mypb.updateDensity();
	mypb.setOperator();

	//for(int ii=0; ii<5; ii++) mypb.solve1step();
	//mypb.solve();
	mypb.solve1step();

	mypb.fpb.printVec(mypb.psi);

	ofstream outFile;
    outFile<<setiosflags(ios::left)<<setiosflags(ios::fixed);
    outFile.open("out.txt");

    int ii;
    for(ii=0; ii<mypb.fpb.Nz; ii++){
    	outFile << mypb.fpb.zs[ii] << "\t" << mypb.psi[ii] << endl;
    }

    outFile.close();
    outFile<<resetiosflags(ios::left)<<resetiosflags(ios::fixed);
	

    return 0;
}