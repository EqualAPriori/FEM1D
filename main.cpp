#include <iostream>
#include "FEM.h"
using namespace std;

int main(){
    FEM f(3);
    cout << "size " << f.ws.size() << " test " << f.xis.at(1) << " weight " << f.ws.at(1) << endl;
    cout << "h " << f.hs[0][0] << ' ' << f.hs[0][1] << ' ' << f.hs[1][0] << ' ' << f.hs[1][1] << endl;


/*
    vector<double> a = {1,2,3,4,5};
    vector<double> b(5,-1);
    cout << b[1] << endl;

    b = a;
    cout << b[1] << endl;
*/

    vector<double> zs = {0,1,2,3,4,5};
    f.setGrid(zs);

    for(int ii = 0; ii <5; ii++){
    	cout << f.zsOfElemE[ii][0] << " " << f.zsOfElemE[ii][1] << " " << f.zsOfElemE[ii][2] << " " << f.zsOfElemE[ii][3] << endl;
    }
    cout << endl;
    for(int ii = 0; ii <5; ii++){
    	cout << f.indOfElemE[ii][0] << " " << f.indOfElemE[ii][1] << " " << f.indOfElemE[ii][2] << " " << f.indOfElemE[ii][3] << endl;
    }


    return 0;
}

