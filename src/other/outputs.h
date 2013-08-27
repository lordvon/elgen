void writeMultiBlockGrid(Grid* grid,char * name){
	//Writes out formatted (ASCII) grid in whole, multi-block PLOT3D format.
	FILE * file = fopen(name, "w");
	//Write total block count.
	fprintf(file,"%d\n",grid->tb);//Number of blocks.
	//Write grid dimensions.
	int b;
	for(b=0;b<grid->tb;b++){
		fprintf(file,"%d\t%d\t%d\n",grid->xdim[b][0],grid->xdim[b][1],grid->xdim[b][2]);
	}
	//Write grids.
	int xi;
	int i,j,k;
	for(b=0;b<grid->tb;b++){
		for (k=0;k<grid->xdim[b][2];k++) {
			for (j=0;j<grid->xdim[b][1];j++){
				for (i=0;i<grid->xdim[b][0];i++){
					xi=i+j*grid->xdim[b][0];
					fprintf(file,"\t%f",grid->x[b][xi]);
				}
				fprintf(file,"\n");
			}
			fprintf(file,"\n");
		}
		for (k=0;k<grid->xdim[b][2];k++) {
			for (j=0;j<grid->xdim[b][1];j++){
				for (i=0;i<grid->xdim[b][0];i++){
					xi=i+j*grid->xdim[b][0];
					fprintf(file,"\t%f",grid->y[b][xi]);
				}
				fprintf(file,"\n");
			}
			fprintf(file,"\n");
		}
		for (k=0;k<grid->xdim[b][2];k++) {
			for (j=0;j<grid->xdim[b][1];j++){
				for (i=0;i<grid->xdim[b][0];i++){
					fprintf(file,"\t%f",0.0);
				}
				fprintf(file,"\n");
			}
			fprintf(file,"\n");
		}
		fprintf(file,"\n");
	}
	fclose(file);
}
void subWriteSoln(FILE *file,double *soln,int b,Grid* g){
	int i,j,k,xi;
	for (k=0;k<g->xdim[b][2];k++){
		for (j=0;j<g->xdim[b][1];j++){
			for (i=0;i<g->xdim[b][0];i++){
				xi=i+j*g->xdim[b][0];
				fprintf(file,"\t%f",soln[xi]);
			}
			fprintf(file,"\n");
		}
		fprintf(file,"\n");
	}
}
void writeMultiBlockGridTranspose(Grid* grid,char * name){
	//Instead of constant-eta lines, constant-xi lines are used.
	//Writes out formatted (ASCII) grid in whole, multi-block PLOT3D format.
	FILE * file = fopen(name, "w");
	//Write total block count.
	fprintf(file,"%d\n",grid->tb);//Number of blocks.
	//Write grid dimensions.
	int b;
	for(b=0;b<grid->tb;b++){
		fprintf(file,"%d\t%d\t%d\n",grid->xdim[b][0],grid->xdim[b][1],grid->xdim[b][2]);
	}
	//Write grids.
	int xi;
	int i,j,k;
	for(b=0;b<grid->tb;b++){
		for (k=0;k<grid->xdim[b][2];k++) {
			for (i=0;i<grid->xdim[b][0];i++){
			for (j=0;j<grid->xdim[b][1];j++){
					xi=i+j*grid->xdim[b][0];
					fprintf(file,"\t%f",grid->x[b][xi]);
				}
				fprintf(file,"\n");
			}
			fprintf(file,"\n");
		}
		for (k=0;k<grid->xdim[b][2];k++) {
			for (i=0;i<grid->xdim[b][0];i++){
			for (j=0;j<grid->xdim[b][1];j++){
					xi=i+j*grid->xdim[b][0];
					fprintf(file,"\t%f",grid->y[b][xi]);
				}
				fprintf(file,"\n");
			}
			fprintf(file,"\n");
		}
		for (k=0;k<grid->xdim[b][2];k++) {
			for (i=0;i<grid->xdim[b][0];i++){
			for (j=0;j<grid->xdim[b][1];j++){
					fprintf(file,"\t%f",0.0);
				}
				fprintf(file,"\n");
			}
			fprintf(file,"\n");
		}
		fprintf(file,"\n");
	}
	fclose(file);
}
void writeMultiBlockCustomSolution(char *name,Grid* g,double **xxsolution){
	//Writes out formatted (ASCII) solution in whole, multi-block PLOT3D format.
	FILE * file = fopen(name, "w");
	int b;
	//Information headers.
	fprintf(file,"%d\n",g->tb);
	for (b=0;b<g->tb;b++){
		fprintf(file,"%d\t%d\t%d\n",g->xdim[b][0],g->xdim[b][1],g->xdim[b][2]);
	}
	//dummy field definitions.
	double **zero=(double **)malloc(sizeof(double)*g->tb);
	for (b=0;b<g->tb;b++){
		zero[b]=(double *)malloc(sizeof(double)*g->xdim[b][3]);
		initv(zero[b],g->xdim[b][3],0);
	}
	//Write fields.
	double mach=0, alpha=0, reyn=0, time=0;
	for (b=0;b<g->tb;b++){
		fprintf(file,"%f\t%f\t%f\t%f\n",mach,alpha,reyn,time);
		subWriteSoln(file,xxsolution[b],b,g);
		subWriteSoln(file,zero[b],b,g);
		subWriteSoln(file,zero[b],b,g);
		subWriteSoln(file,zero[b],b,g);
		subWriteSoln(file,zero[b],b,g);
	}
	fclose(file);
}
