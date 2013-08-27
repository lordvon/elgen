int blockMatch(int xi1,int xi2,int b1,int b2,
		Grid* g){
	double matchTol=1e-6;
	double diff;
	int match=0;
	diff=	sqrt(
				pow(g->x[b1][xi1]-g->x[b2][xi2],2)+
				pow(g->y[b1][xi1]-g->y[b2][xi2],2)
			);
	if(diff<matchTol){ match=1; }
	return match;
}
int getSideXdim(int side, int *xdim){
	int sidedim;
	if((side==0) | (side==2)){
		sidedim=xdim[0];
	} else {
		sidedim=xdim[1];
	}
	return sidedim;
}
int getInterfaceBlock(Grid* g,int block,int side){
	//from block and side, deduce adjacent block.
	int otherblock,matchedblock;
	int xi1,xi2;
	if(side==0){ xi1=0; }
	else if(side==1){ xi1=g->xdim[block][0]-1; }
	else if(side==2){ xi1=g->xdim[block][3]-g->xdim[block][0]; }
	else if(side==3){ xi1=0; }
	matchedblock=-1;
	//search through all other blocks to find interfaced block.
	for(otherblock=0;otherblock<g->tb;otherblock++){
		if(side==0){ xi2=g->xdim[otherblock][3]-g->xdim[otherblock][0]; }
		else if(side==1){ xi2=0; }
		else if(side==2){ xi2=0; }
		else if(side==3){ xi2=g->xdim[otherblock][0]-1; }
		if(blockMatch(xi1,xi2,block,otherblock,g)>0){
			matchedblock=otherblock;
			break;
		}
	}
	if(0){
		printf("Error: Interfacing block not found at for block %d (side %d)!\n",block,side);
	} else if(0){//matchedblock>-1) {
		printf("Matched block for block %d (side %d): %d\n",block,side,matchedblock);
	}
	return matchedblock;
}
int getInterfaceSide(Grid*g,int b1,int b2){
	int matched=-1;
	int s;
	for(s=0;s<4;s++){
		if(b2==getInterfaceBlock(g,b1,s)){
			matched=s;
			break;
		}
	}
	if(matched<0){
		printf("Error: Interfacing side not found between blocks %d and %d!\n",b1,b2);
	}
	return matched;
}


int getTransverseDim(Grid*g,int nb,int*oldblocks,int directionIndex){
	int b,oldb,total;
	total=0;
	for(b=0;b<nb;b++){
		oldb=oldblocks[b];
		total+=g->xdim[oldb][directionIndex];
	}
	total-=(nb-1);
	return total;
}
void copyMultiBlockGrid02(Grid*g,Grid*gm,int side,int nb,int*oldblocks,int newb){
	int i,j,xi,base;
	int oldb;
	int jmax,jmin,xd0,xd1,xd3;
	int currentBlock;

	//Strategy for 0,2 joins.
	base=0;
	currentBlock=0;
	printf("%d, %d %d\n",nb,oldblocks[0],oldblocks[1]);
	while(currentBlock<nb){
		oldb=oldblocks[currentBlock];
		xd0=g->xdim[oldb][0];
		xd1=g->xdim[oldb][1];
		xd3=g->xdim[oldb][3];
		//printf("old block for newb %d: %d\n",newb,oldb);
		jmin=0;
		jmax=xd1-1;
		if(currentBlock<nb-1){
			switch(side){
			case 0: jmin++; break;
			case 2: jmax--; break;
			}
		}
		for(j=jmin;j<=jmax;j++){
			for(i=0;i<=xd0-1;i++){
				xi=i+j*xd0;
				//printf("%g, %g\n",g->x[oldb][xi],g->y[oldb][xi]);
				gm->x[newb][base+xi]=g->x[oldb][xi];
				gm->y[newb][base+xi]=g->y[oldb][xi];
				//printf("%g, %g\n",gm->x[newb][base+xi],gm->y[newb][base+xi]);
			}
		}
		base+=xd3-xd0;
		currentBlock++;
	}
}
void copyMultiBlockGrid13(Grid*g,Grid*gm,int side,int nb,int*oldblocks,int newb){
	int b,i,j,xi1,xi2;
	int oldb;
	int imax,imin;
	int xd0,xd1;
	int currentBlock;
	int *istart=malloc(sizeof(int)*nb);
	//make istart.
	int total=0;
	for(b=0;b<nb;b++){
		oldb=oldblocks[b];
		istart[b]=total;
		total+=(g->xdim[oldb][0]-1);
	}
	//Get new dimensions.
	//Strategy for 1,3 joins.
	currentBlock=0;
	//printf("%d %d %d\n",istart[0],istart[1],istart[2]);
	//printf("%d %d %d\n",oldblocks[0],oldblocks[1],oldblocks[2]);
	while(currentBlock<nb){
		oldb=oldblocks[currentBlock];
		xd0=g->xdim[oldb][0];
		xd1=g->xdim[oldb][1];
		imin=0;
		imax=g->xdim[oldb][0]-1;
		if(currentBlock<nb-1){
			switch(side){
			case 1: imax--; break;
			case 3: imin++; break;
			}
		}
		for(j=0;j<=xd1-1;j++){
			for(i=imin;i<=imax;i++){
				xi1=i+j*xd0;
				xi2=istart[currentBlock]+i+j*gm->xdim[newb][0];
				gm->x[newb][xi2]=g->x[oldb][xi1];
				gm->y[newb][xi2]=g->y[oldb][xi1];
			}
		}
		currentBlock++;
	}
	free(istart);
}
void copyMultiBlockGrid(Grid*g,Grid*gm,int side,int nb,int*oldblocks,int newb){
	if(side%2==0){ copyMultiBlockGrid02(g,gm,side,nb,oldblocks,newb); }
	else { copyMultiBlockGrid13(g,gm,side,nb,oldblocks,newb); }
}
void mergedGrid(Grid*g,Grid*gm,Merges*m){
	int ntb,ntn,mid,side,bid0,bid1,sidedim;
	int bc,bb;
	int oldb,oldb0,newb;
	int *oldblocks,nob;
	//Initialize mask.
	m->mask=malloc(sizeof(int)*g->tb);
	for(oldb=0;oldb<g->tb;oldb++){ m->mask[oldb]=-1; }
	//fill mask, new grid total dimensions, merge sides.
	ntb=g->tb;
	ntn=g->tn;
	for(mid=0;mid<m->tm;mid++){
		nob=m->nb[mid];
		oldblocks=m->merges[mid];
		bid0=oldblocks[0];
		bid1=oldblocks[1];
		side=getInterfaceSide(g,bid0,bid1);
		ntb-=(nob-1);
		if(side%2==0){ sidedim=g->xdim[bid0][0]; }
		else { sidedim=g->xdim[bid0][1]; }
		ntn-=(nob-1)*sidedim;
		m->sides[mid]=side;
		for(oldb=0;oldb<nob;oldb++){ m->mask[oldblocks[oldb]]=mid; }
	}
	gm->tb=ntb;
	gm->tn=ntn;
	//malloc gm fields.
	gm->xdim=mbinfo(gm->tb);
	gm->x=malloc(sizeof(double*)*gm->tb);
	gm->y=malloc(sizeof(double*)*gm->tb);
	m->key=malloc(sizeof(int*)*gm->tb);
	m->keynb=malloc(sizeof(int)*gm->tb);
	//fill key.
	int mask,newxd3;
	bc=0;
	int *checklist=malloc(sizeof(int)*m->tm);
	for(mid=0;mid<m->tm;mid++){ checklist[mid]=0; }
	for(oldb=0;oldb<g->tb;oldb++){
		mask=m->mask[oldb];
		if(mask<0){
			m->keynb[bc]=1;
			m->key[bc]=malloc(sizeof(int));
			m->key[bc][0]=oldb;
			gm->xdim[bc][0]=g->xdim[oldb][0];
			gm->xdim[bc][1]=g->xdim[oldb][1];
			gm->xdim[bc][2]=g->xdim[oldb][2];
			gm->xdim[bc][3]=g->xdim[oldb][3];
			bc++;
		} else {
			if(checklist[mask]!=1){
				nob=m->nb[mask];
				oldblocks=m->merges[mask];
				m->keynb[bc]=nob;
				m->key[bc]=malloc(sizeof(int)*nob);
				newxd3=0;
				for(bb=0;bb<nob;bb++){
					m->key[bc][bb]=oldblocks[bb];
					newxd3+=g->xdim[m->key[bc][bb]][3];
				}
				side=m->sides[mask];
				if(side%2==0){
					sidedim=g->xdim[m->key[bc][0]][0];
					gm->xdim[bc][0]=sidedim;
					gm->xdim[bc][1]=getTransverseDim(g,nob,oldblocks,1);
				} else {
					sidedim=g->xdim[m->key[bc][0]][1];
					gm->xdim[bc][0]=getTransverseDim(g,nob,oldblocks,0);
					gm->xdim[bc][1]=sidedim;
				}
				newxd3-=(nob-1)*sidedim;
				gm->xdim[bc][2]=1;
				gm->xdim[bc][3]=newxd3;
				gm->x[bc]=malloc(sizeof(double)*newxd3);
				gm->y[bc]=malloc(sizeof(double)*newxd3);
				checklist[mask]=1;
				bc++;
			}
		}
	}
	//fill new grid.
	for(mid=0;mid<m->tm;mid++){ checklist[mid]=0; }
	for(newb=0;newb<gm->tb;newb++){
		oldb0=m->key[newb][0];
		mid=m->mask[oldb0];
		if(mid<0){
			gm->x[newb]=g->x[oldb0];
			gm->y[newb]=g->y[oldb0];
		} else if(checklist[mid]<1){
			//printf("newb: %d\n",newb);
			side=m->sides[mid];
			copyMultiBlockGrid(g,gm,side,m->nb[mid],m->merges[mid],newb);
			checklist[mid]=1;
		}
	}
	free(checklist);
}
void extractBlockVert(Grid*g,Grid*gm,int oldbi,int*oldblocks,int newb){
	int bcounter,jstart,xd0;
	int i,j,xistart,xi1,xi2;

	int oldb=oldblocks[oldbi];

	xd0=gm->xdim[newb][0];
	jstart=0;
	if(oldbi>0){
		for(bcounter=0;bcounter<oldbi;bcounter++){
			jstart+=g->xdim[oldblocks[bcounter]][1]-1;
		}
		jstart++;
	}
	xistart=jstart*xd0;

	for(j=0;j<g->xdim[oldb][1];j++){
		for(i=0;i<xd0;i++){
			xi1=i+j*xd0;
			xi2=xistart+xi1;
			g->x[oldb][xi1]=gm->x[newb][xi2];
			g->y[oldb][xi1]=gm->y[newb][xi2];
		}
	}
}
void extractBlockHori(Grid*g,Grid*gm,int oldbi,int*oldblocks,int newb){
	int bcounter,istart,xd1;
	int i,j,xi1,xi2;
	int oldxd0,newxd0;

	int oldb=oldblocks[oldbi];

	xd1=gm->xdim[newb][1];
	istart=0;
	if(oldbi>0){
		for(bcounter=0;bcounter<oldbi;bcounter++){
			istart+=g->xdim[oldblocks[bcounter]][0]-1;
		}
		//istart++;
	}

	oldxd0=g->xdim[oldb][0];
	newxd0=gm->xdim[newb][0];
	for(j=0;j<xd1;j++){
		for(i=0;i<oldxd0;i++){
			xi1=i+j*oldxd0;
			xi2=istart+i+j*newxd0;
			g->x[oldb][xi1]=gm->x[newb][xi2];
			g->y[oldb][xi1]=gm->y[newb][xi2];
		}
	}
}
void extractBlock(Grid*g,Grid*gm,int side,int oldbi,int*oldblocks,int newb){
	if(side==2){
		extractBlockVert(g,gm,oldbi,oldblocks,newb);
	} else if(side==1){
		extractBlockHori(g,gm,oldbi,oldblocks,newb);
	} else {
		printf("Error: sides 0 and 3 are not valid merging locations.\n");
	}
}
void splitGrid(Grid*g,Grid*gm,Merges*m){
	int newb,oldb,bcounter,side,mid;
	for(newb=0;newb<gm->tb;newb++){
		oldb=m->key[newb][0];
		if(m->keynb[newb]==1){
			g->x[oldb]=gm->x[newb];
			g->y[oldb]=gm->y[newb];
		} else {
			mid=m->mask[oldb];
			side=m->sides[mid];
			for(bcounter=0;bcounter<m->keynb[newb];bcounter++){
				extractBlock(g,gm,side,bcounter,m->key[newb],newb);
			}
		}
	}
}
