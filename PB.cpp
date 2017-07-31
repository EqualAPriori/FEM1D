#include "PB.h"


PB::PB():fpb(2){};
PB::PB(int ord){
	order = ord; 
	fpb.setOrder(ord); 
}

//default, create equally spaced grid
void PB::setGrid(double Len, int Nx){
	L = Len;
	Nz = Nx;
	dx = Len/((double)Nx-1.0);

	nodes.resize(Nz);
	int ii;
	for(ii=0; ii<Nz; ii++){
		nodes[ii] = dx * (double) ii;
	}

	fpb.setGrid(nodes);
	//fpb.printVec(fpb.zs);
	//fpb.printVec2D(fpb.zsOfElemE);

	myDGESV.refresh(fpb.Nz);
}
void PB::setGrid(vector<double>& coords){
	nodes.assign(coords.begin(),coords.end());
	fpb.setGrid(nodes);
	//fpb.printVec2D(fpb.zsOfElemE);
	myDGESV.refresh(fpb.Nz);
}

void PB::bulkInit(){
	int iSpecies, iz; 
	psi.resize(fpb.Nz);
	fill(psi.begin(),psi.end(),0.0);

	kap20 = 0.0; 
	for(iSpecies=0; iSpecies<Nspecies; iSpecies++){
		kap20+=rho0[iSpecies]*q[iSpecies]*q[iSpecies];
	}
	kap20*=4*M_PI*lb;

	rho.resize(Nspecies, vector<double>(fpb.Nz));

	netchg.resize(fpb.Nz);
	fill(netchg.begin(),netchg.end(),0.0);
	kap2.resize(fpb.Nz);
	fill(kap2.begin(),kap2.end(),0.0);
	netchgLessDH.resize(fpb.Nz);
	fill(netchgLessDH.begin(),netchgLessDH.end(),0.0);
	//for(auto &v: rho){fill(rho0.begin(),rho0.end(),0.0);}
	//for(auto &v: rho){memset(&v[0], 0, sizeof(v[0]) * v.size());}

	for(iSpecies=0; iSpecies<Nspecies; iSpecies++){
    	fill(rho[iSpecies].begin(), rho[iSpecies].end(), rho0[iSpecies]);
    	for(iz=0; iz<fpb.Nz; iz++){
    		kap2[iz] += 4*M_PI*lb*rho[iSpecies][iz]*q[iSpecies]*q[iSpecies];
    		netchg[iz] += q[iSpecies]*rho[iSpecies][iz];
    	}
	}

	errflag = 0; //since psi = 0, can't calculate error (denominator is zero)
	err = 1.0;
	iterations = 0;
}

void PB::updateDensity(){
	int iSpecies, iz;
	fill(kap2.begin(),kap2.end(),0.0);
	fill(netchg.begin(),netchg.end(),0.0);
	for(iSpecies=0; iSpecies<Nspecies; iSpecies++){
		for(iz=0; iz<fpb.Nz; iz++){
			rho[iSpecies][iz] = rho0[iSpecies] * exp( -q[iSpecies]*psi[iz] );
			kap2[iz] += 4*M_PI*lb*rho[iSpecies][iz]*q[iSpecies]*q[iSpecies];
			netchg[iz] += q[iSpecies]*rho[iSpecies][iz];
		}
	}
	for(iz=0; iz<fpb.Nz; iz++){
    	netchgLessDH[iz] = 4*M_PI*lb*netchg[iz] + kap20*psi[iz];
   	}
}

void PB::setOperator(){
	fpb.BuildLaplacian();
	//fpb.printSqMat(fpb.D2);

	fpb.BuildMass(kap20); 	//default is linearized local density approximation
	//fpb.BuildMass(kap2);
	//fpb.printVec(fpb.M);

	fpb.combineOperator();
	//fpb.printSqMat(fpb.Op);

	fpb.BuildBoundaryOperator(BCs);
	//cout << fpb.BCs[0][0] << " " << fpb.BCs[0][1] << " " << fpb.BCs[0][2] << endl;
}

void PB::setBC(double const bc[2][3]){
	BCs[0][0] = bc[0][0];
	BCs[0][1] = bc[0][1];
	BCs[0][2] = bc[0][2];
	BCs[1][0] = bc[1][0];
	BCs[1][1] = bc[1][1];
	BCs[1][2] = bc[1][2];
}

// 1 Step of the iteration
void PB::solve1step(){
	fpb.BuildSource(netchgLessDH);
	//fpb.BuildSource(0.0);
	//fpb.printVec(fpb.f);

	//setOperator();  //shoot, LAPACK modifies my operator matrix, have to recreate it each time!!!
	//fpb.BuildBoundary(BCs);
	fpb.BuildBoundarySource();
	//fpb.printVec(fpb.f);
	//fpb.printVec(fpb.f);
	//fpb.printSqMat(fpb.Op);

	//Calculate solution using LAPACK -- use my wrapper function instead
	double *A = &fpb.Op[0];
	double *b = &fpb.f[0];
	vector<double> psiNew; psiNew.resize(fpb.Nz,0.0);
	myDGESV.solve(A,b,&psiNew[0]);

	//calculate error and get next iteration
	if(errflag==1){
		err = calcError(psiNew,psi);
	}else{
		errflag = 1;
		err = 1.0;
	}
	cout << "Calculation Error: " << err << endl;

	nextSol(psiNew);
	updateDensity();
}

// solve until less than tolerance.
void PB::solve(){
	while(iterations < iterMax && err > tol){
		solve1step();
		iterations++;
	}
	if(iterations >= iterMax){
		cout << "Failed convergence, iterations EXCEEDED iterMax, but error is: " << err << endl;
	}else if(err < tol){
		cout << "Successful convergence, " << iterations << " iterations: with error: " << err << endl;
	}

}


double PB::calcError(vector<double>& ynew, vector<double>& yold){
	if(yold.size() != ynew.size()){
		cout << "Error: new and old vectors of incommensurate size" << endl;
		throw;
	}
	int nn = yold.size();
	int iz;
	double num = 0.0;
	double denom = 0.0;;
	for(iz=0; iz<nn; iz++){
		num += (ynew[iz]-yold[iz])*(ynew[iz]-yold[iz]);
		denom += yold[iz]*yold[iz];
	}
	return sqrt(num/denom);
}

void PB::nextSol(vector<double>& psiNew, int flag){ //calculates next psi and stores it
	if(flag==0){ //Picard. mixing rate determined by error
		int nn = psi.size();
		int iz;
		//mix = 0.1;
		//mix = mixBase/(mixBase + err);
		for(iz=0; iz<nn; iz++){
			psi[iz] = (1-mix)*psi[iz] + mix* (psiNew)[iz];
		}
	}else{
		cout << "Error: unknown iteration method. use (0:Picard), (1:DIIS)" << endl;
		throw;
	}
	//fpb.printVec(psiNew);
}

