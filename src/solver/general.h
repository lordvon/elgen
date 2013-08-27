void mallocSolver(Grid*g,Solver*s){
	int b,xd3;
	s->dx=malloc(sizeof(double*)*g->tb);
	s->dy=malloc(sizeof(double*)*g->tb);
	for(b=0;b<g->tb;b++){
		xd3=g->xdim[b][3];
		s->dx[b]=malloc(sizeof(double)*xd3);
		s->dy[b]=malloc(sizeof(double)*xd3);
	}
}
double computeRMSResidual(Grid*g,Solver*s){
	//After the computation the residuals matrix can be reset.
	int b,xd3,xi;
	double runningTotal=0;
	for(b=0;b<g->tb;b++){
		xd3=g->xdim[b][3];
		for(xi=0;xi<xd3;xi++){
			runningTotal+=(pow(s->dx[b][xi],2)+pow(s->dy[b][xi],2));
		}
	}
	runningTotal/=g->tn;
	double residual=sqrt(runningTotal);
	return residual;
}
