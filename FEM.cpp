// Class file for FEM
#include "FEM.h"
#include <iostream>

FEM::FEM(int ord){
	if(ord<1){
		cout << "order is faulty, defaulting to linear order" << endl;
		order = 1;
	}else if(ord>5){
		cout << "order exceeds current implementation, defaulting to max order = 5" << endl;
		order = 5;
	}else{
		order = ord;
	}
	NGLL = ord+1;

	GetGLL(NGLL, & xis, & ws, & hs);
}

// Establish FEM grid given coarse nodes (i.e. depending on GLL order, FEM fills in intermediate points)
// zcoords should be in strictly increasing order, elementary validation only done in xiToZ() and zToXi()
void FEM::setGrid(vector<double>& zcoords){
	zcoarse.assign(zcoords.begin(), zcoords.end());
	Nel = zcoarse.size() - 1;
	Nz = Nel*order + 1;

	zs.resize(Nz,0);
	zsOfElemE.resize(Nel, vector<double>(NGLL,0));
	indOfElemE.resize(Nel, vector<int>(NGLL,0));

	int el, ii;
	vector<double> temp;
	for(el = 0; el < Nel; ++el){
		//temp = xiToZ(xis,zcoarse[el],zcoarse[el+1]);
		xiToZ(xis,zcoarse[el],zcoarse[el+1], &temp);
		zsOfElemE[el].assign(temp.begin(),temp.end());

		for(ii=0; ii<NGLL; ++ii){
			indOfElemE[el][ii] = (el)*(NGLL-1) + ii;
		}
	}
}

// Gauss-Lobatto-Legendre quadrature values. Taken from Ampuero (Caltech) and Wolfram
// http://mathworld.wolfram.com/LobattoQuadrature.html
// 2017.07.25
void FEM::GetGLL(int npts, vector<double>* xi, vector<double>* w, vector<vector<double> >* h){
	xi->resize(npts,0);
	w->resize(npts,0);
	h->resize(npts, vector<double>(npts,0));

	if(npts < 2){
		cout << "order is out of range, defaulting to linear order" << endl;
		npts = 2;
	}else if(npts > 6){
		cout << "order exceeds current implementation, defaulting to max order=5 (6pts)" << endl;
		npts = 6;
	}

	if(npts == 2){
		(*xi) = {-1.0, 1.0};
		(*w) = {1.0, 1.0};
		(*h) = { {-0.5, -0.5},
				 { 0.5,  0.5} };
	}else if(npts == 3){
		(*xi) = {-1.0, 0.0, 1.0};
		(*w) = {1.0/3.0, 4.0/3.0, 1.0/3.0};
		(*h) = {{-1.5, -0.5, 0.5},
				{2.0, 0.0, -2.0},
				{-0.5, 0.5, 1.5}};
	}else if(npts == 4){
		(*xi) = {-1.0, -sqrt(5.0)/5.0, sqrt(5.0)/5.0, 1.0};
		(*w) = {1.0/6.0, 5.0/6.0, 5.0/6.0, 1.0/6.0};
		(*h) = {{-3.000000000000000, -0.8090169943749475, 0.3090169943749474, -0.5000000000000000},
 				{4.045084971874738, 0.000000000000000E+00, -1.118033988749895, 1.545084971874737},
 				{-1.545084971874737, 1.118033988749895, 0.000000000000000E+00, -4.045084971874738},
 				{0.5000000000000000, -0.3090169943749474, 0.8090169943749475, 3.000000000000000}};
	}else if(npts == 5){
		(*xi) = {-1.0, -sqrt(21.0)/7.0, 0.0, sqrt(21.0)/7.0, 1.0};
		(*w) = {0.1, 49.0/90.0, 32.0/45.0, 49.0/90.0, 0.1};
		(*h) = {{-5.000000000000000, -1.240990253030983, 0.3750000000000000, -0.2590097469690172, 0.5000000000000000},
 				{6.756502488724241, 0.000000000000000E+00, -1.336584577695453, 0.7637626158259734, -1.410164177942427},
 				{-2.666666666666667, 1.745743121887939, 0.000000000000000E+00, -1.745743121887940, 2.666666666666667},
 				{1.410164177942427, -0.7637626158259734, 1.336584577695453, 0.000000000000000E+00, -6.756502488724238},
 				{-0.5000000000000000, 0.2590097469690171, -0.3750000000000000, 1.240990253030982, 5.000000000000000}};
	}else if(npts == 6){
		(*xi) = {-1.0, sqrt( (7.0 + 2.0*sqrt(7.0))/21.0 ), sqrt( (7.0 - 2.0*sqrt(7.0))/21.0 ), \
					   sqrt( (7.0 - 2.0*sqrt(7.0))/21.0 ), sqrt( (7.0 + 2.0*sqrt(7.0))/21.0 ), 1.0};
		(*w) = {1.0/15.0, (14.0-sqrt(7.0))/30.0, (14.0+sqrt(7.0))/30.0, \
						  (14.0+sqrt(7.0))/30.0, (14.0-sqrt(7.0))/30.0, 1.0/15.0};	
		(*h) = {{ -7.500000000000000, -1.786364948339096, 0.4849510478535692, -0.2697006108320389, 0.2377811779842314, -0.5000000000000000},
 				{10.14141593631967, 0.000000000000000E+00, -1.721256952830233, 0.7863566722232407, -0.6535475074298002, 1.349913314190488},
 				{-4.036187270305349, 2.523426777429455, 0.000000000000000E+00, -1.752961966367866, 1.152828158535930, -2.244684648176167},
 				{2.244684648176167, -1.152828158535930, 1.752961966367866, 0.000000000000000E+00, -2.523426777429455, 4.036187270305349},
 				{-1.349913314190488, 0.6535475074298002, -0.7863566722232407, 1.721256952830233, 0.000000000000000E+00, -10.14141593631967},
 				{0.5000000000000000, -0.2377811779842314, 0.2697006108320389, -0.4849510478535692, 1.786364948339096, 7.500000000000000}};	
	}
}


//converting coordinates
void FEM::xiToZ(vector<double>& xi, double zleft, double zright, vector<double>* z){
		z->assign(xi.begin(),xi.end());
		if(zleft >= zright){
			cout << "Error: zleft >= zright in xiToZ()" << endl;
			return;
		}

		double len = fabs(zright - zleft);
		int ii;
		for(ii=0; ii < xi.size(); ii++){
			(*z)[ii] = zleft + len*(xi[ii] + 1.0)/2.0;
		}
	}

void FEM::zToXi(vector<double>& z, double zleft, double zright, vector<double>* xi){
		xi->assign(z.begin(),z.end());
		if(zleft >= zright){
			cout << "Error: zleft >= zright in zToXi()" << endl;
			return;
		}

		double len = fabs(zright - zleft);
		int ii;
		for(ii=0; ii < z.size(); ii++){
			(*xi)[ii] = 2.0*(z[ii] - zleft)/len - 1.0;
		}
	}


//int newData[5] = {0, 1, 2, 3, 4};

// assign(#,val) function will reinitialize the vector to # elements with value `val`
// assign(newData, newData + count); where it will range over the array from first pointer to last

// for 2D MxN:
//std::vector<std::vector<int>> fog(
//    M,
//    std::vector<int>(N, valToFillMatrices));

//matrix.resize( row_count , vector<int>( column_count , initialization_value ) );

//myVector.resize(n);
//for (int i = 0; i < n; ++i)
//    myVector[i].resize(m);

