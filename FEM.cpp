// Class file for FEM class
//
// Note currently coded s.t. a FEM class instance can be built to build matrices/write results 
// to arrays whos addresses are pointed in (instead of opearating only on its own matrices).
// "Normal" class behavior (operating on own data) is provided via overloaded constructors
// although if feeding in user matrices, one needs to be careful that they are properly zeroed first
// since I only implement sparse zeroing.
//

#include "FEM.h"

FEM::FEM(int ord){
	setOrder(ord);
}
void FEM::setOrder(int ord){
	if(ord<1){
		cout << "order is too low, defaulting to linear order" << endl;
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
// 2017.07.25
void FEM::setGrid(vector<double>& zcoords){
	int flag = 0; //to see if dim of zcoord changed
	if(zcoarse.size() != zcoords.size())
		flag=1;

	zcoarse.assign(zcoords.begin(), zcoords.end());
	Nel = zcoarse.size() - 1;
	Nz = Nel*order + 1;

	if(flag==1){ //currently misses edge case when `order` changes. hopefully never change order in usage
		zs.resize(Nz,0);
		zsOfElemE.resize(Nel, vector<double>(NGLL));
		indOfElemE.resize(Nel, vector<int>(NGLL));
	}

	int el, ii, ind;
	vector<double> temp;
	for(el = 0; el < Nel; ++el){
		//temp = xiToZ(xis,zcoarse[el],zcoarse[el+1]);
		xiToZ(xis,zcoarse[el],zcoarse[el+1], &temp);
		zsOfElemE[el].assign(temp.begin(),temp.end());

		for(ii=0; ii<NGLL; ++ii){
			ind = (el)*(NGLL-1) + ii;
			indOfElemE[el][ii] = ind;
			zs[ind] = zsOfElemE[el][ii];
		}
	}

	if(flag==1){
		D2loc.resize(NGLL, vector<double>(NGLL));
		D2.resize(Nz*Nz,0); fill(D2.begin(),D2.end(),0.0);
		M.resize(Nz*1,0);	// by construction a diagonal matrix, so we only keep the diagonal.
		K.resize(Nz*Nz,0);	fill(K.begin(),K.end(),0.0);
		Op.resize(Nz*Nz,0); fill(Op.begin(),Op.end(),0.0);
		f.resize(Nz,0);
	}
}

// Building Laplacian Matrix
// default is overloaded function with eps = 1 everywhere
// gives NEGATIVE laplacian -d(eps (du)), or -(eps*u'(z))' like in Poisson equation
// 2017.07.26
void FEM::BuildLaplacian(vector<double>& eps, vector<double>* D2out, bool clear){
	if(eps.size() != zs.size()){
		cout << "Error: `epsilon` size is incongruent with zs" << endl;
		throw;
	}

	int el,ii,jj,kk,ind1,ind2,indColMajor;
	double jacob;

	//sparse zero-resetting. must be done separately from assembly becuase there are repeats
	//other elements were already zeroed in setGrid();
	//reset ONLY if boolean flag=1 (default)
	if(clear==1){
		for(el=0; el<Nel; el++){
			for(ii=0; ii<NGLL; ii++){
				ind1 = indOfElemE[el][ii];
				for(jj=0; jj<NGLL; jj++){
					ind2 = indOfElemE[el][jj];
					(*D2out)[ij(ind1,ind2,Nz)] = 0;
				}
			}	
		}	
	}	

	//assembly step
	for(el=0; el<Nel; el++){
		jacob = 2.0/(zcoarse[el+1]-zcoarse[el]);

		for(ii=0; ii<NGLL; ii++){
			ind1 = indOfElemE[el][ii];
			for(jj=0; jj<NGLL; jj++){
				ind2 = indOfElemE[el][jj];
				indColMajor = ij(ind1,ind2,Nz);
				for(kk=0; kk<NGLL; kk++){ //index over nodes in element
					(*D2out)[indColMajor] += jacob * hs[ii][kk] * hs[jj][kk] * ws[kk] * eps[indOfElemE[el][kk]];
				}
			}
		}
	}
}
void FEM::BuildLaplacian(vector<double>* D2out, bool clear, double eps ){
	int el,ii,jj,kk,ind1,ind2,indColMajor;
	double jacob;

	//cout << "clear flag: " << clear << " scalar eps: " << eps << endl;

	//sparse zero-resetting. must be done separately from assembly becuase there are repeats
	//other elements were already zeroed in setGrid();
	//reset ONLY if boolean flag=1 (default)
	if(clear==1){
		for(el=0; el<Nel; el++){
			for(ii=0; ii<NGLL; ii++){
				ind1 = indOfElemE[el][ii];
				for(jj=0; jj<NGLL; jj++){
					ind2 = indOfElemE[el][jj];
					(*D2out)[ij(ind1,ind2,Nz)] = 0;
				}
			}
		}
	}

	//assembly step
	for(el=0; el<Nel; el++){
		jacob = 2.0/(zcoarse[el+1]-zcoarse[el]);

		for(ii=0; ii<NGLL; ii++){
			ind1 = indOfElemE[el][ii];
			for(jj=0; jj<NGLL; jj++){
				ind2 = indOfElemE[el][jj];
				indColMajor = ij(ind1,ind2,Nz);
				for(kk=0; kk<NGLL; kk++){ //index over nodes in element
					(*D2out)[indColMajor] += jacob * hs[ii][kk] * hs[jj][kk] * ws[kk] * eps;
				}
			}
		}
	}
}


// Building FEM source
// 2017.07.27
void FEM::BuildSource(vector<double>& source, vector<double>* fout){
	if(source.size() != zs.size()){
		cout << "Error: `source` size is incongruent with zs" << endl;
		throw;
	}

	int el,ii,indZ;
	double jacob;

	//Unfortunately need to reset vector
    //&memset(&f[0], 0, f.size() * sizeof f[0]); // fastest in all circumstances
	fill(fout->begin(),fout->end(),0); //must use -O3 optimization flag for speed! see https://stackoverflow.com/questions/8848575/fastest-way-to-reset-every-value-of-stdvectorint-to-0

	// Assemble source vector
	for(el=0; el<Nel; el++){
		for(ii=0; ii<NGLL; ii++){
			indZ = indOfElemE[el][ii];
			jacob = 0.5*(zcoarse[el+1]-zcoarse[el]);
			(*fout)[indZ] += source[indZ] * ws[ii] *jacob;
		}
	}
}
void FEM::BuildSource(double source, vector<double>* fout){
	fout->resize(Nz);
	fill(fout->begin(),fout->end(),0);

	int el,ii,indZ;
	double jacob;

	// Assemble source vector
	for(el=0; el<Nel; el++){
		for(ii=0; ii<NGLL; ii++){
			indZ = indOfElemE[el][ii];
			jacob = 0.5*(zcoarse[el+1]-zcoarse[el]);
			(*fout)[indZ] += source * ws[ii] *jacob;
		}
	}
}


// Building FEM mass (should be diagonal term if using GLL evaluation
// 2017.07.27
void FEM::BuildMass(vector<double>& mass, vector<double>* Mout){
	if(mass.size() != zs.size()){
		cout << "Error: `source` size is incongruent with zs" << endl;
		throw;
	}

	int el, ii, indZ;
	double jacob;
	// Only a sparse reset is needed, along the diagonal
	for(ii=0; ii<Nz; ii++){
		(*Mout)[ii] = 0;
	}

	// Assemble mass diagonal
	for(el=0; el<Nel; el++){
		jacob = 0.5*(zcoarse[el+1]-zcoarse[el]);
		for(ii=0; ii<NGLL; ii++){
			indZ = indOfElemE[el][ii];

			//indColMajor = ij(indZ,indZ,Nz);
			//(*Mout)[indColMajor] += mass[indZ] * ws[ii] * jacob;
			(*Mout)[indZ] += mass[indZ] * ws[ii] * jacob;
		}
	}
}
void FEM::BuildMass(double mass, vector<double>* Mout){
	Mout->resize(Nz);
	fill(Mout->begin(),Mout->end(),0);
	int el, ii, indZ;
	double jacob;

	for(el=0; el<Nel; el++){
		jacob = 0.5*(zcoarse[el+1]-zcoarse[el]);
		for(ii=0; ii<NGLL; ii++){
			indZ = indOfElemE[el][ii];

			//indColMajor = ij(indZ,indZ,Nz);
			//(*Mout)[indColMajor] += mass[indZ] * ws[ii] * jacob;
			(*Mout)[indZ] += mass * ws[ii] * jacob;
		}
	}
}
//builds into NxN matrix, *WITHOUT* clearing
// 2017.07.31
void FEM::BuildMassMatrix(vector<double>& mass, vector<double>* OpOut){
	if(mass.size() != zs.size()){
		cout << "Error: `source` size is incongruent with zs" << endl;
		throw;
	}

	int el, ii, indZ, indColMajor;
	double jacob;

	// Assemble mass onto diagonal of output matrix
	for(el=0; el<Nel; el++){
		jacob = 0.5*(zcoarse[el+1]-zcoarse[el]);
		for(ii=0; ii<NGLL; ii++){
			indZ = indOfElemE[el][ii];
			indColMajor = ij(indZ,indZ,Nz);
			(*OpOut)[indColMajor] += mass[indZ] * ws[ii] * jacob;
		}
	}
} 
void FEM::BuildMassMatrix(double mass, vector<double> *OpOut){
	int el, ii, indZ, indColMajor;
	double jacob;

	// Assemble mass onto diagonal of output matrix
	for(el=0; el<Nel; el++){
		jacob = 0.5*(zcoarse[el+1]-zcoarse[el]);
		for(ii=0; ii<NGLL; ii++){
			indZ = indOfElemE[el][ii];
			indColMajor = ij(indZ,indZ,Nz);
			(*OpOut)[indColMajor] += mass * ws[ii] * jacob;
		}
	}
}

// Building FEM boundary, modifies both the operator matrix `A` and the source term `b` for Ax=b
// done in a way s.t. for Dirichlet, keep extra rows and columns, don't delete row/column
// In the future, can consider allowing for internal interfaces/nodal conditions
// Format of boundary matrix (2x3) is for: eps*dpsi/dx - alpha*psi = beta (note sign in front of alpha)
// <position>	<type>		<alpha> 	<beta>
// left bdry 	0 = dirich.    #		   #
// right bdry   1 = mixed 	0 if Neumann   #
//							1 default if dirich.
// ASSUMES `A` is square, and `b` is the same size
// 2017.07.27
void FEM::BuildBoundary(double const bc[2][3], vector<double>* A, vector<double>* b){
	BCs[0][0] = bc[0][0];
	BCs[0][1] = bc[0][1];
	BCs[0][2] = bc[0][2];
	BCs[1][0] = bc[1][0];
	BCs[1][1] = bc[1][1];
	BCs[1][2] = bc[1][2];

	int icol;
	int N = sqrt(A->size());

	if(N!=b->size()){
		cout << "Matrix `A` size incommensurate with forcing term `b`" << endl;
		throw;
	}
	//Left
	if(bc[0][0]==0){//Dirichlet
		(*b)[0] = bc[0][2] * (zcoarse[1]-zcoarse[0]); //scale by length of element s.t. commensurate with scale of F.E.M.
		(*A)[0] = 1.0 * (zcoarse[1]-zcoarse[0]);		 //this scaling only makes sense if on same zgrid as FEM...
		for(icol=1; icol< N; icol++){
			(*A)[ij(0,icol,N)] = 0.0;
		}
	}else{//Mixed
		(*b)[0] -= bc[0][2]; // -beta
		(*A)[ij(0,0,N)] += bc[0][1]; // +alpha
	}
	//Right
	if(bc[1][0]==0){//Dirichlet
		(*b)[N-1] = bc[1][2] * (zcoarse[Nel]-zcoarse[Nel-1]);
		(*A)[ij(N-1,N-1,N)] = 1.0 * (zcoarse[Nel]-zcoarse[Nel-1]);
		for(icol=0; icol<N-1; icol++){
			(*A)[ij(N-1,icol,N)] = 0.0;
		}
	}else{//Mixed
		(*b)[N-1] += bc[1][2]; // +beta
		(*A)[ij(N-1,N-1,N)] -= bc[1][1]; // -alpha
	}	
}

//Build only the operator to reflect BCs
void FEM::BuildBoundaryOperator(double const bc[2][3], vector<double>* A){
	BCs[0][0] = bc[0][0];
	BCs[0][1] = bc[0][1];
	BCs[0][2] = bc[0][2];
	BCs[1][0] = bc[1][0];
	BCs[1][1] = bc[1][1];
	BCs[1][2] = bc[1][2];

	int icol;
	int N = sqrt(A->size());

	//Left
	if(bc[0][0]==0){//Dirichlet
		(*A)[0] = 1.0 * (zcoarse[1]-zcoarse[0]);		 //this scaling only makes sense if on same zgrid as FEM...
		for(icol=1; icol< N; icol++){
			(*A)[ij(0,icol,N)] = 0.0;
		}
	}else{//Mixed
		(*A)[ij(0,0,N)] += bc[0][1];
	}
	//Right
	if(bc[1][0]==0){//Dirichlet
		(*A)[ij(N-1,N-1,N)] = 1.0 * (zcoarse[Nel]-zcoarse[Nel-1]);
		for(icol=0; icol<N-1; icol++){
			(*A)[ij(N-1,icol,N)] = 0.0;
		}
	}else{//Mixed
		(*A)[ij(N-1,N-1,N)] -= bc[1][1];
	}		
}

//Modify the source term to reflect BCs, using previously calculated BC
//this way avoids repeatedly modifying the operator (which is wrong)
void FEM::BuildBoundarySource(vector<double>* b){
	int N = sqrt(Op.size());
	if(b->size() != N){
		cout << "Error: Source term's size incommensurate with operator Op stored in FEM object." << endl;
		throw;
	}

	if(BCs[0][0]==0){//Dirichlet
		(*b)[0] = BCs[0][2] * (zcoarse[1]-zcoarse[0]); //scale by length of element s.t. commensurate with scale of F.E.M.
	}else{//Mixed
		(*b)[0] -= BCs[0][2];
	}
	//Right
	if(BCs[1][0]==0){//Dirichlet
		(*b)[N-1] = BCs[1][2] * (zcoarse[Nel]-zcoarse[Nel-1]);
	}else{//Mixed
		(*b)[N-1] += BCs[1][2];
	}	
}

// Building FEM integral
// need to feed in a matrix of the kernel *evaluated* at the nodes
// ASSUMES input matrices are square and on same grid as stored in FEM object
// zeroes the `Kout` matrix first.
// 2017.07.27
void FEM::BuildIntegral(vector<double>& kernel, vector<double> *Kout){
	if(kernel.size() != Kout->size()){
		cout << "Error: `kernel` size is incongruent with `Kout`" << endl;
		throw;
	}
	if(kernel.size() != Nz*Nz){
		cout << "Error: `kernel` size incongruent with `Nz`" << endl;
	}
	if(Kout->size() != Nz*Nz){
		cout << "Error: `Kout` size incongruent with `Nz`" << endl;
	}

	int el1,el2, ii,jj, ind1,ind2, indColMajor;
	double w1,w2, l1,l2;

	// zeroing
	for(ii=0; ii<Kout->size(); ii++){ (*Kout)[ii] = 0; }

	//assembly step
	for(el1=0; el1<Nel; el1++){
		for(el2=0; el2<Nel; el2++){
			for(ii=0; ii<NGLL; ii++){
				for(jj=0; jj<NGLL; jj++){
					ind1 = indOfElemE[el1][ii];
					ind2 = indOfElemE[el2][jj];

					w1 = ws[ii];
					w2 = ws[jj];
					l1 = zcoarse[el1+1] - zcoarse[el1];
					l2 = zcoarse[el2+1] - zcoarse[el2];

					indColMajor = ij(ind1,ind2,Nz);

					//cout << ind1 << " " << ind2 << " " << w1 << " " << w2 << " " << l1 << " " << l2 << endl << flush;
					(*Kout)[indColMajor] += 0.25*w1*w2*l1*l2*kernel[indColMajor];
				}
			}
		}
	}
}

// Gather pieces into final operator `A`
// 2017.07.27
void FEM::combineOperator(vector<double>& Lap, vector<double>& mass, vector<double>* A){
	A->assign(Lap.begin(),Lap.end());
	int ii,indColMajor;
	int Npts = mass.size();
	for(ii=0; ii<Npts; ii++){
		indColMajor = ij(ii,ii,Npts);
		(*A)[indColMajor] += mass[ii];
	}
}

// Building Integro-Differential operator directly into an operator, 
// to eliminate unnecessary clearing and copying
// 2017.07.31
void FEM::BuildIntDiffOperator(vector<double> &kernel, \
	vector<double> &eps, vector<double> &mass, vector<double> *OpOut){

	int clear = 0;
	BuildIntegral(kernel, OpOut);
	BuildLaplacian(eps, OpOut, clear);
	BuildMassMatrix(mass, OpOut);
}//vector eps, vector mass
void FEM::BuildIntDiffOperator(vector<double> &kernel, \
	vector<double> &eps, double mass, vector<double> *OpOut){

	int clear = 0;
	BuildIntegral(kernel, OpOut);
	BuildLaplacian(eps, OpOut, clear);
	BuildMassMatrix(mass, OpOut);
}//vector eps, scalar mass
void FEM::BuildIntDiffOperator(vector<double> &kernel, \
	double eps, vector<double> &mass, vector<double> *OpOut){

	int clear = 0;
	BuildIntegral(kernel, OpOut);
	BuildLaplacian(OpOut, clear, eps); //unfortunately my Laplacian function with scalar eps has eps as last argument... (typically=1.0, can ignore)
	BuildMassMatrix(mass, OpOut);
}//scalar eps, vector mass
void FEM::BuildIntDiffOperator(vector<double> &kernel, \
	double eps, double mass, vector<double> *OpOut){

	int clear = 0;
	BuildIntegral(kernel, OpOut);
	BuildLaplacian(OpOut, clear, eps); //unfortunately my Laplacian function with scalar eps has eps as last argument... (typically=1.0, can ignore)
	BuildMassMatrix(mass, OpOut);
}//scalar eps, scalar mass


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

// converting coordinates
// 2017.07.25
void FEM::xiToZ(vector<double>& xi, double zleft, double zright, vector<double>* z){
		if(z->size() != xi.size())
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
		if(xi->size() != z.size())
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


