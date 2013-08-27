void phisub(int j,int b,Grid*g,Transformation*t,Equation*e){
	int xi,i;
	int xd0=g->xdim[b][0];
	double diff;
	for(i=0;i<xd0;i++){
		xi=i+j*xd0;
		diff=fabs(t->xxi[b][xi])-fabs(t->yxi[b][xi]);
		if(diff>0){
			e->phi[b][xi]=-t->xxixi[b][xi]/t->xxi[b][xi];
		} else {
			e->phi[b][xi]=-t->yxixi[b][xi]/t->yxi[b][xi];
		}
	}
}
void psisub(int i,int b,Grid*g,Transformation*t,Equation*e){
	int xi,j;
	int xd0=g->xdim[b][0];
	int xd1=g->xdim[b][1];
	double diff;
	for (j=0;j<xd1;j++){
		xi=i+j*xd0;
		diff=fabs(t->xeta[b][xi])-fabs(t->yeta[b][xi]);
		if(diff>0){
			e->psi[b][xi]=-t->xetaeta[b][xi]/t->xeta[b][xi];
		} else {
			e->psi[b][xi]=-t->yetaeta[b][xi]/t->yeta[b][xi];
		}
	}
}
void phisubcorner(int j,int b,Grid*g,Transformation*t,Equation*e,Numerics*n){
	//assumes phi is filled already.
	int xi,i;
	int xd0=g->xdim[b][0];
	int check1,check2;
	double sign;
	int iv=1;
	for(i=1;i<xd0-1;i++){
		xi=i+j*xd0;
		check1=fabs(e->phi[b][xi+iv]/e->phi[b][xi])<n->cornertol;
		check2=fabs(e->phi[b][xi-iv]/e->phi[b][xi])<n->cornertol;
		if(check1 & check2){
			printf("Phi spike limited at block %d (%d, %d)\n",b,i,j);
			sign=e->phi[b][xi]/fabs(e->phi[b][xi]);
			e->phi[b][xi]=sign*n->spikelim;
		}
	}
}
void psisubcorner(int i,int b,Grid*g,Transformation*t,Equation*e,Numerics*n){
	//assumes psi is filled already.
	int xi,j;
	int xd0=g->xdim[b][0];
	int xd1=g->xdim[b][1];
	int check1,check2;
	double sign;
	int iv=xd0;
	for (j=1;j<xd1-1;j++){
		xi=i+j*xd0;
		check1=fabs(e->psi[b][xi+iv]/e->psi[b][xi])<n->cornertol;
		check2=fabs(e->psi[b][xi-iv]/e->psi[b][xi])<n->cornertol;
		if(check1 & check2){
			printf("Psi spike limited at block %d (%d, %d)\n",b,i,j);
			sign=e->psi[b][xi]/fabs(e->psi[b][xi]);
			e->psi[b][xi]=sign*n->spikelim;
		}
	}
}
void phipsi(Grid*g,Transformation*t,Equation*e,Numerics*n){
	//Computes control terms.
	int b,i,j,xi;
	int xd0,xd1,xd3;
	int xi1,xi2;
	double diff;
	for(b=0;b<g->tb;b++){
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		xd3=g->xdim[b][3];
		phisub(0		,b,g,t,e);
		phisubcorner(0		,b,g,t,e,n);
		phisub(xd1-1	,b,g,t,e);
		phisubcorner(xd1-1	,b,g,t,e,n);
		for (i=0;i<xd0;i++){
			xi1=i;
			xi2=xi1+xd3-xd0;
			diff=e->phi[b][xi2]-e->phi[b][xi1];
			for(j=1;j<xd1-1;j++){
				xi=i+j*xd0;
				e->phi[b][xi]=e->phi[b][xi1]+j*diff/(xd1-1);
			}
		}
		psisub(0		,b,g,t,e);
		psisubcorner(0		,b,g,t,e,n);
		psisub(xd0-1	,b,g,t,e);
		psisubcorner(xd0-1	,b,g,t,e,n);
		for (j=0;j<xd1;j++){
			xi1=j*xd0;
			xi2=xi1+xd0-1;
			diff=e->psi[b][xi2]-e->psi[b][xi1];
			for(i=1;i<xd0-1;i++){
				xi=i+j*xd0;
				e->psi[b][xi]=e->psi[b][xi1]+i*diff/(xd0-1);
			}
		}
	}
}
