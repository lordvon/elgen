void inputs(Grid*g,Merges*m,Numerics*n){
	readNumerics(n,"in/numerics");
	prepGrid(g,"in/blockDict");
	prepMerge(m,"in/mergeDict");
}

void initialization(Grid*g,Transformation*t,Equation*e,Solver*s,Numerics*n){
	initTransformation(t,g);
	mallocEquation(g,e);
	phipsi(g,t,e,n);
	mallocSolver(g,s);
}
void update(Grid*g,Transformation*t,Equation*e){
	inverseMetrics(g,t);
	fillCoefficients(g,t,e);
}
void solve(Grid*g,Equation*e,Solver*s){
	tridiagonalSolverThomasSerial(g,e,g->x,e->c,e->d,s->dx);
	tridiagonalSolverThomasSerial(g,e,g->y,e->gamma,e->delta,s->dy);
}
