#include "Green.h"


Green::Green():fe(2){};
Green::Green(int ord){
	order = ord; 
	fe.setOrder(ord); 
}

//default, create equally spaced grid
void Green::setGrid(double Len, int Nx){
	L = Len;
	Nz = Nx;
	dx = Len/((double)Nx-1.0);

	nodes.resize(Nz);
	int ii;
	for(ii=0; ii<Nz; ii++){
		nodes[ii] = dx * (double) ii;
	}

	fe.setGrid(nodes);
	//fe.printVec(fe.zs);
	//fe.printVec2D(fe.zsOfElemE);

	myDGESV.refresh(fe.Nz);

	//Green object variables
	psi.resize(fe.Nz);
	diagKernelPrev.resize(fe.Nz);
}
void Green::setGrid(vector<double>& coords){
	nodes.assign(coords.begin(),coords.end());
	fe.setGrid(nodes);
	//fe.printVec2D(fe.zsOfElemE);
	myDGESV.refresh(fe.Nz);

	//Green object variables
	psi.resize(fe.Nz);
	diagKernelPrev.resize(fe.Nz);
	fill(diagKernelPrev.begin(), diagKernelPrev.end(), 0.0);
}

void Green::setBC(double const bc[2][3]){
	BCs[0][0] = bc[0][0];
	BCs[0][1] = bc[0][1];
	BCs[0][2] = bc[0][2];
	BCs[1][0] = bc[1][0];
	BCs[1][1] = bc[1][1];
	BCs[1][2] = bc[1][2];
}
// right now default is for 1 wall value only 
// (even if symmetric walls, solve on half space and still only left bdry modified)
void Green::setlb(double lbs, double lbw){
	lb = lbs;
	lbwall = lbw;
	BCs[0][1] = lbw/lbs; //leaving right bdry either Dirichlet or Homogeneous Neumann
}

// Building the Bare Operator (without the diagonal wavenumber k-term)
// diagonal flag is whether the kernel is pt-charge delta functions (i.e. diagonal)
// may also be helpful to have an additional ftn that handles just delta-type kernel contributions
// operator is stored in the FEM object
//
// note that, for delta functions, using Kernel matrix gives slightly different result ]
// than Kernel diagonal vector. Because when using Kernel matrix, if one is not careful to scale 
// the delta function at interface correctly to account for the fact that is only HALF of a 'tent', 
// then one will end up missing ~half of the value of the delta function.
// 
// in other words, the boundary elements of the Kernel *VECTOR* is NOT just the diagonal of the
// Kernel *MATRIX*
//
// for equally spaced terms, the K(0,0) term of the Kernel Matrix will have a 0.25 factor,
// whereas the K(0)=K(0,0) term of the Kernel diagonal Vector will have a 0.5 factor 
// --> off by factor of 2!
//
void Green::BuildBareOperator(vector<double> &kernel, double eps){
	bool diagonal;
	if(kernel.size() == fe.Nz){
		diagonal = 1; //kernel is just delta ftns
	}else if(kernel.size() == fe.Nz*fe.Nz){
		diagonal = 0;
	}else{
		cout << "Error: kernel has inappropriate size." << endl;
		throw;
	}

	if(diagonal==0){//full operator
		fe.BuildIntDiffOperator(kernel, eps, 0.0);
	}else{//sparse (in fact, diagonal) integral kernel
		int clear=0;
		fill(fe.Op.begin(), fe.Op.end(), 0.0); //zero out the operator on `Op`
		fe.BuildLaplacian(&fe.Op, clear, eps); //build Laplacian on top of `Op`
		modifyDeltaFtnKernel(kernel, &fe.Op); //build Integration kernel on top of `Op`
	}

	fe.BuildBoundaryOperator(BCs);
}
void Green::BuildBareOperator(vector<double> &kernel, vector<double> &eps){
	bool diagonal;
	if(kernel.size() == fe.Nz){
		diagonal = 1; //kernel is just delta ftns
	}else if(kernel.size() == fe.Nz*fe.Nz){
		diagonal = 0;
	}else{
		cout << "Error: kernel has inappropriate size." << endl;
		throw;
	}

	if(diagonal==0){//full operator
		fe.BuildIntDiffOperator(kernel, eps, 0.0);
	}else{//sparse (in fact, diagonal) integral kernel
		int clear=0; //
		fill(fe.Op.begin(), fe.Op.end(), 0.0); //zero out the operator on `Op`
		fe.BuildLaplacian(eps, &fe.Op, clear); //build Laplacian on top of `Op`
		modifyDeltaFtnKernel(kernel, &fe.Op); //build Integration kernel on top of `Op`
	}

	fe.BuildBoundaryOperator(BCs);
};

//note this ftn does NOT clear the diagonal.
// written this way, kernel is just the value of the scalar multiplying delta function
// NOT in the form with a value ~1/dx/w factor of an approximate delta function.
// this function will take care of things automatically.
void Green::modifyDeltaFtnKernel(vector<double> &kernel, vector<double> *Op){
	fe.BuildMassMatrix(diagKernelPrev, Op); //first undo the previous kernel.
	fe.BuildMassMatrix(kernel, Op);
	int ii = 0;
	for(ii=0; ii<kernel.size(); ii++){
		diagKernelPrev[ii] = -kernel[ii]; //store the NEGATIVE of the kernel for undoing later on.
	}
}
void Green::modifyDeltaFtnKernel( double kernel, vector<double> *Op){
	fe.BuildMassMatrix(diagKernelPrev, Op);
	fe.BuildMassMatrix(kernel, Op);
	int ii=0;
	fill(diagKernelPrev.begin(), diagKernelPrev.end(), -kernel);
}

// Modifying just the wave-number contribution (which is diagonal!)
// object keeps track of previous k used s.t. only have to modify diagonal instead of completely
// rebuilding the entire dense, square matrix
void Green::kOperator(double k, vector<double> *Op){
	fe.BuildMassMatrix(kprev, Op);
	fe.BuildMassMatrix(k, Op);
	kprev = - k; //store the NEGATIVE of the wavenumber `k` for undoing later on.
}


// Linear Solve, given source delta `z0` and wavenumber `k`
// since my positions `zs` are sorted, can locate index with log(Nz) time (binary search).
void Green::solve(int loc, double k){
	kOperator(k); //
	//Construct the delta function representation
	buildDelta(loc);

	fe.BuildBoundarySource();
	//fe.printVec(fe.f);
	//fe.printVec(fe.f);
	//fe.printSqMat(fe.Op);

	//Calculate solution using LAPACK -- use my wrapper function instead
	double *A = &fe.Op[0];
	double *b = &fe.f[0];
	//vector<double> psiNew; psiNew.resize(fe.Nz,0.0);
	myDGESV.solve(A,b,&psi[0]);
}
void Green::solve(double z0, double k){
	int loc = find(z0);
	solve(loc, k);
}


// Finding the magnitude of the delta function located at a node.
// as argument takes INDEX into the desired node.
void Green::buildDelta(int loc, vector<double> *f){
	fill((*f).begin(), (*f).end(), 0.0);
	(*f)[loc] = 1.0;
}
// version where cartesian COORDINATE of node is given. use binary search to find index fast.
// alternative is to use a recursive binarySearch function
void Green::buildDelta(double z0, vector<double> *f){
	if( z0 < fe.zs[0] || z0>fe.zs[fe.zs.size()-1]){
		cout << "Error: delta source location `z0` out of bounds" << endl;
	}

	int loc = find(z0,fe.zs);

	//cout << z0 << " closest to " << fe.zs[loc] << " at index " << loc << endl;
	buildDelta(loc, f);
}


// search for index of coordinate `z0` in `pos`
double Green::find(double z0, vector<double> &pos){
	int Nz = pos.size();
	int lowbd = 0;
	int upbd = Nz;
	int loc;

	// a log(N) index search
	loc = Nz/2;
	while(upbd - lowbd > 1){
		if( z0 == pos[loc]){
			lowbd = loc;
			upbd = loc;
			break;
		}else if( z0 > pos[loc]){
			lowbd = loc;
		}else{
			upbd = loc;
		}
		loc = lowbd + (upbd-lowbd)/2;
	}

	if(upbd-lowbd == 1){
		//cout << "z0 " << z0 << " chosen between grid points, ROUNDING (0.5 --> 1.0)" << endl;
		if( (z0-pos[lowbd]) >= (pos[upbd]-z0) ){
			//cout << "rounding up" << endl;
			loc = upbd;
		}else{
			//cout << "rounding down" << endl;
			loc = lowbd;
		}
	}
	cout << "using position " << pos[loc] << " (input was " << z0 << ")" << endl;
	return loc;
}


