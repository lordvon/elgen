void countTotalBlocks(Grid* g,char *inputfile){
	int maxchars=200;
	int blocknumber,totalblocks;
	char line[maxchars];
	char filename[maxchars];
	FILE * file = fopen (inputfile, "rt");
	int success;
	//Read total number of blocks.
	totalblocks=0;
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %s",&blocknumber,filename);
		if(line[0]=='#'){ break; }
		if(success>0){ totalblocks++; }
	}
	g->tb=totalblocks;
	fclose(file);
	g->blockNames=charmalloc(g->tb,200);
}
void countTotalMerges(Merges* ms,char *inputfile){
	int maxchars=200;
	int totalmerges;
	int mid,bs,bcounter,b;
	char line[maxchars];
	FILE * file = fopen (inputfile, "rt");
	int success;
	totalmerges=0;
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %d",&mid,&bs);
		if(line[0]=='#'){ break; }
		if(success>0){
			totalmerges++;
			bcounter=0;
			while(fgets(line, maxchars, file) != NULL){
				success = sscanf (line, "%d",&b);
				if(success>0){
					bcounter++;
					if(bcounter>=bs){ break; }
				}
			}
		}
	}
	ms->tm=totalmerges;
	ms->merges=malloc(sizeof(int*)*ms->tm);
	ms->nb=malloc(sizeof(int)*ms->tm);
	ms->sides=malloc(sizeof(int)*ms->tm);
	fclose(file);
}
void readMergeDict(Merges* ms,char *inputfile){
	int maxchars=200;
	int mid,bs,b,bcounter;
	char line[maxchars];
	FILE * file = fopen (inputfile, "rt");
	int success;
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %d",&mid,&bs);
		if(line[0]=='#'){ break; }
		if(success>0){
			ms->nb[mid]=bs;
			ms->merges[mid]=malloc(sizeof(int)*bs);
			bcounter=0;
			while(fgets(line, maxchars, file) != NULL){
				success = sscanf (line, "%d",&b);
				if(success>0){
					ms->merges[mid][bcounter]=b;
					bcounter++;
					if(bcounter>=bs){ break; }
				}
			}
		}
	}
	fclose(file);
}
void readBlockDict(Grid* g,char *inputfile){
	int maxchars=200;
	int blocknumber;
	char line[maxchars];
	char filename[maxchars];
	FILE * file = fopen (inputfile, "rt");
	int success;
	//Read file names for each block.
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %s",&blocknumber,filename);
		if(line[0]=='#'){ break; }
		if(success>0){
			strcpy(g->blockNames[blocknumber],filename);
		}
	}
	fclose(file);
}
void readGrids(Grid* g){
	g->xdim=mbinfo(g->tb);
	g->x=(double **)malloc((g->tb)*sizeof(double *));
	g->y=(double **)malloc((g->tb)*sizeof(double *));
	int b,success,dim[3],maxchars=200;
	char line[maxchars];
	double xcoord,ycoord;
	g->tn=0;
	for(b=0;b<g->tb;b++){
		//readSingleBlockGrid(grid->blockNames[b],grid->xdim[b],grid->x[b],grid->y[b]);
		FILE * file = fopen (g->blockNames[b], "rt");
		//Get grid dimensions.
		while(fgets(line, maxchars, file) != NULL)
		{
			success = sscanf (line, "%d %d %d", &dim[0], &dim[1], &dim[2]);
			if(success>0){
				//printf("%s dimensions: %d %d %d.\n",grid->blockNames[b],dim[0],dim[1],dim[2]);
				break;
			}
		}
		g->xdim[b][0]=dim[0];
		g->xdim[b][1]=dim[1];
		g->xdim[b][2]=dim[2];
		g->xdim[b][3]=g->xdim[b][0]*g->xdim[b][1]*g->xdim[b][2];
		g->tn+=g->xdim[b][3];
		//Read grid.
		g->x[b]=(double *)malloc(sizeof(double)*g->xdim[b][3]);
		g->y[b]=(double *)malloc(sizeof(double)*g->xdim[b][3]);
		int xi=0;
		while(fgets(line, maxchars, file) != NULL) {
			success = sscanf (line, "%lf %lf", &xcoord, &ycoord);
			if(success>0){
				if(xi>=g->xdim[b][3]){ printf("Warning: More grid points were available than specified.\n"); }
				g->x[b][xi]=xcoord;
				g->y[b][xi]=ycoord;
				xi++;
			}
		}
		fclose(file);
	}
}
void prepGrid(Grid* g,char *blockFileName){
	countTotalBlocks(g,blockFileName);
	readBlockDict(g,blockFileName);
	readGrids(g);
}
void prepMerge(Merges*m,char*mergeFileName){
	countTotalMerges(m,mergeFileName);
	readMergeDict(m,mergeFileName);
}
