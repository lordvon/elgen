
void mallocTransformation(Transformation*t,Grid*g){
	int b,xd3;
	t->xxi=malloc(sizeof(double*)*g->tb);
	t->xeta=malloc(sizeof(double*)*g->tb);
	t->yxi=malloc(sizeof(double*)*g->tb);
	t->yeta=malloc(sizeof(double*)*g->tb);
	t->xxixi=malloc(sizeof(double*)*g->tb);
	t->xetaeta=malloc(sizeof(double*)*g->tb);
	t->yxixi=malloc(sizeof(double*)*g->tb);
	t->yetaeta=malloc(sizeof(double*)*g->tb);
	for(b=0;b<g->tb;b++){
		xd3=g->xdim[b][3];
		t->xxi[b]=malloc(sizeof(double)*xd3);
		t->xeta[b]=malloc(sizeof(double)*xd3);
		t->yxi[b]=malloc(sizeof(double)*xd3);
		t->yeta[b]=malloc(sizeof(double)*xd3);
		t->xxixi[b]=malloc(sizeof(double)*xd3);
		t->xetaeta[b]=malloc(sizeof(double)*xd3);
		t->yxixi[b]=malloc(sizeof(double)*xd3);
		t->yetaeta[b]=malloc(sizeof(double)*xd3);
	}
}
void inverseMetricsSub(int xi,int iv,double sign,
		double**metx1,double**mety1,double**metx2,double**mety2,
		Grid*g,int b){
	metx1[b][xi]=sign*0.5*(3*g->x[b][xi]-4*g->x[b][xi+iv]+g->x[b][xi+2*iv]);
	mety1[b][xi]=sign*0.5*(3*g->y[b][xi]-4*g->y[b][xi+iv]+g->y[b][xi+2*iv]);
	metx2[b][xi]=sign*(g->x[b][xi]-2*g->x[b][xi+iv]+g->x[b][xi+2*iv]);
	mety2[b][xi]=sign*(g->y[b][xi]-2*g->y[b][xi+iv]+g->y[b][xi+2*iv]);
}
void inverseMetrics(Grid*g,Transformation*t){
	//Computes 1st and 2nd order FD derivatives of x and y wrt xi and eta
	//second-order accurate, except for one-sided 2nd order derivatives.
	int b,i,j,xi;
	int xd0,xd1;
	int i1,i2,j1,j2;
	for(b=0;b<g->tb;b++){
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		//interior
		for (j=0;j<xd1;j++){
			for (i=1;i<xd0-1;i++){
				xi=i+j*xd0;
				t->xxi[b][xi]=0.5*(g->x[b][xi+1]-g->x[b][xi-1]);
				t->yxi[b][xi]=0.5*(g->y[b][xi+1]-g->y[b][xi-1]);
				t->xxixi[b][xi]=g->x[b][xi+1]-2*g->x[b][xi]+g->x[b][xi-1];
				t->yxixi[b][xi]=g->y[b][xi+1]-2*g->y[b][xi]+g->y[b][xi-1];
			}
		}
		for (j=1;j<xd1-1;j++){
			for (i=0;i<xd0;i++){
				xi=i+j*xd0;
				t->xeta[b][xi]=0.5*(g->x[b][xi+xd0]-g->x[b][xi-xd0]);
				t->yeta[b][xi]=0.5*(g->y[b][xi+xd0]-g->y[b][xi-xd0]);
				t->xetaeta[b][xi]=g->x[b][xi+xd0]-2*g->x[b][xi]+g->x[b][xi-xd0];
				t->yetaeta[b][xi]=g->y[b][xi+xd0]-2*g->y[b][xi]+g->y[b][xi-xd0];
			}
		}
		//boundaries
		i1=0;i2=xd0-1;
		for (j=0;j<xd1;j++){
			inverseMetricsSub(i1+j*xd0,1,-1,
					t->xxi,t->yxi,t->xxixi,t->yxixi,g,b);
			inverseMetricsSub(i2+j*xd0,-1,1,
					t->xxi,t->yxi,t->xxixi,t->yxixi,g,b);
		}
		j1=0;j2=xd1-1;
		for (i=0;i<xd0;i++){
			inverseMetricsSub(i+j1*xd0,xd0,-1,
					t->xeta,t->yeta,t->xetaeta,t->yetaeta,g,b);
			inverseMetricsSub(i+j2*xd0,-xd0,1,
					t->xeta,t->yeta,t->xetaeta,t->yetaeta,g,b);
		}
	}
}
void initTransformation(Transformation*t,Grid*g){
	mallocTransformation(t,g);
	inverseMetrics(g,t);
}
