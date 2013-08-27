

void fillA(Grid*g,Transformation*t,Equation*e){
	int b,xi,xd3;
	for(b=0;b<g->tb;b++){
		xd3=g->xdim[b][3];
		for (xi=0;xi<xd3;xi++){
			e->A1[b][xi]=pow(t->xeta[b][xi],2)+pow(t->yeta[b][xi],2);
			e->A2[b][xi]=t->xxi[b][xi]*t->xeta[b][xi]+
					t->yxi[b][xi]*t->yeta[b][xi];
			e->A3[b][xi]=pow(t->xxi[b][xi],2)+pow(t->yxi[b][xi],2);
		}
	}
}
void fillab(Grid*g,Equation*e){
	//Computes coefficients of the i-1 and i+1 terms in the Thompson equation (a/alpha, b/beta in the project description).
	int b,j,xi;
	int xd0,xd1,xd3;
	for(b=0;b<g->tb;b++){
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		xd3=g->xdim[b][3];
		for(xi=0;xi<xd3;xi++){
			e->a[b][xi]=e->A1[b][xi]*(1+0.5*e->phi[b][xi]);
			e->b[b][xi]=e->A1[b][xi]*(1-0.5*e->phi[b][xi]);
		}
		//Boundary Condition
		for (j=0;j<xd1;j++){
			xi=j*xd0;
			e->a[b][xi]=0;
			e->b[b][xi]=0;
			xi=j*xd0+xd0-1;
			e->a[b][xi]=0;
			e->b[b][xi]=0;
		}
	}
}
void fillc(Grid*g,Equation*e,double**c,double**grid){
	//Computes coefficients that dont depend on positions at j at time n+1 in the Thompson equation (c/gamma in the project description).
	int b,j,xi,xd0,xd1,xd3;
	double term1,term2;
	for(b=0;b<g->tb;b++){
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		xd3=g->xdim[b][3];
		for (xi=0;xi<xd3;xi++){
			term1=0.5*e->A2[b][xi]*(
					grid[b][xi+1+xd0]+
					grid[b][xi-1-xd0]-
					grid[b][xi+1-xd0]-
					grid[b][xi-1+xd0]
					        );
			term2=-e->A3[b][xi]*(
					grid[b][xi+xd0]+
					grid[b][xi-xd0]+
					0.5*e->psi[b][xi]*(	grid[b][xi+xd0]-
										grid[b][xi-xd0])
							);
			c[b][xi]=term1+term2;
		}
		//Boundary Condition
		for (j=0;j<xd1;j++){
			xi=j*xd0;
			c[b][xi]=grid[b][xi];
			xi=j*xd0+xd0-1;
			c[b][xi]=grid[b][xi];
		}
	}
}
void filld(Grid*g,Equation*e){
	//Computes coefficients at i,j term in the Thompson equation (d/delta in the project description).
	int b,j,xi,xd0,xd1,xd3;
	for(b=0;b<g->tb;b++){
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		xd3=g->xdim[b][3];
		for(xi=0;xi<xd3;xi++){
			e->d[b][xi]=-2*(e->A1[b][xi]+e->A3[b][xi]);
		}
		//Boundary Condition
		for (j=0;j<xd1;j++){
			xi=j*xd0;
			e->d[b][xi]=1;
			xi=j*xd0+xd0-1;
			e->d[b][xi]=1;
		}
	}
}
void fillCoefficients(Grid*g,Transformation*t,Equation*e){
	fillA(g,t,e);
	fillab(g,e);
	fillc(g,e,e->c,g->x);
	fillc(g,e,e->gamma,g->y);
	filld(g,e);
	mbcopy(g,e->d,e->delta);
}
