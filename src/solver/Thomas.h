/*
void tridiagonalSolverThomasSerial(double ** grid, double ** a, double ** b, double ** c, double ** d, int imax, int jmax, double ** diff){
	//Iterates once through tridiagonal system for every j-line, iterating j-2 tridiagonal systems in total per function call.
	//Adds to the residuals matrix the square of the difference between the new and old grid.
	int i,j;
	double transient;
	for (j=2;j<jmax;j++){
		//Forward substitution.
		for(i=2;i<=imax;i++){
			transient=b[i][j]/d[i-1][j];
			d[i][j]-=transient*a[i-1][j];
			c[i][j]-=transient*c[i-1][j];
		}
		//Back substitution.
		c[imax][j]/=d[imax][j];
		for(i=imax-1;i>=1;i--){
			c[i][j]=(c[i][j]-a[i][j]*c[i+1][j])/d[i][j];
		}
		for(i=2;i<imax;i++){
			diff[i][j]=c[i][j]-grid[i][j];
			grid[i][j]=c[i][j];
		}
	}
}
*/
void tridiagonalSolverThomasSerial(Grid*g,Equation*e,double**grid,
		double**c,double**d,double ** diff){
	//Iterates once through tridiagonal system for every j-line, iterating j-2 tridiagonal systems in total per function call.
	//Adds to the residuals matrix the square of the difference between the new and old grid.
	int b,i,j,xi,xd0,xd1;
	double transient;
	for(b=0;b<g->tb;b++){
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		for (j=1;j<xd1-1;j++){
			//Forward substitution.
			for(i=1;i<xd0;i++){
				xi=i+j*xd0;
				transient=e->b[b][xi]/d[b][xi-1];
				d[b][xi]-=transient*e->a[b][xi-1];
				c[b][xi]-=transient*c[b][xi-1];
			}
			//Back substitution.
			xi=j*xd0+(xd0-1);
			c[b][xi]/=d[b][xi];
			for(i=xd0-2;i>=0;i--){
				xi=i+j*xd0;
				c[b][xi]=(c[b][xi]-e->a[b][xi]*c[b][xi+1])/d[b][xi];
			}
			for(i=1;i<xd0-1;i++){
				xi=i+j*xd0;
				diff[b][xi]=c[b][xi]-grid[b][xi];
				grid[b][xi]=c[b][xi];
			}
		}
	}
}
