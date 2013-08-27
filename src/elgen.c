#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct {//Grid
	int tb,tn;//totalblocks,totalnodes
	int **xdim;//[b][idim,jdim,kdim,blockNodes]
	double **x,**y;//[b][xi]
	char **blockNames;//for reading in the grid.
} Grid;
typedef struct {//Merges
	int tm;//totalmerges
	int *nb;//index by merge id,number of blocks in each merge.
	int **merges;//[m]{all of bs}
	int *sides;//indexed by merge id; side at which first block in each merge is merged.
	int *mask;//indexed by old block id; >=0: merge id, <0: not in a merge.
	int **key;//indexed by new merged block id -> lists all of the old block ids that are contained.
	int *keynb;//indexed by new merged block id -> number of old blocks in each new block.
} Merges;
typedef struct {//Transformation
	double **xxi,**xeta,**yxi,**yeta;
	double **xxixi,**xetaeta,**yxixi,**yetaeta;
} Transformation;
typedef struct {//Equation
	double **phi,**psi;
	double **A1,**A2,**A3;
	double **a,**b,**c,**d,**gamma,**delta;
	double **xphi,**yphi,**xpsi,**ypsi;
} Equation;
typedef struct {//Solver
	double **dx,**dy;
} Solver;
typedef struct {//Numerics
	int useControlTerms;
	int writeResiduals;
	double tol;
	double difftol;
	double cornertol;
	double spikelim;
} Numerics;

#include "other/dataStructures.h"
#include "other/inputs.h"
#include "other/outputs.h"

#include "computationalDomain/grid.h"
#include "computationalDomain/transformation.h"
#include "computationalDomain/connectivity.h"

#include "equation/general.h"
#include "equation/controlTerms.h"
#include "equation/coefficients.h"

#include "solver/general.h"
#include "solver/Thomas.h"
#include "solver/steps.h"


int main(void) {
	Grid g;
	Grid gm;//merged version
	Merges m;
	Transformation t;
	Equation e;
	Solver s;
	Numerics n;

	//pre
	inputs(&g,&m,&n);
	mergedGrid(&g,&gm,&m);
	initialization(&gm,&t,&e,&s,&n);
	//solve
	int iter=0;double res=1;
	while (res>n.tol){
		iter++;
		update(&gm,&t,&e);
		solve(&gm,&e,&s);
		res=computeRMSResidual(&gm,&s);
		//printf("Iteration %d Residual: %g\n",iter,res);
	}
	printf("Done (%d iterations, final RMS residual: %e).\n",iter,res);

	//post
	splitGrid(&g,&gm,&m);
	writeMultiBlockGrid(&g,"out/g.xyz");
	writeMultiBlockGrid(&gm,"out/gm.xyz");
	writeMultiBlockGridTranspose(&gm,"out/gmcon.xyz");

	/*
	writeMultiBlockCustomSolution("out/phi.q",&gm,e.phi);
	writeMultiBlockCustomSolution("out/psi.q",&gm,e.psi);
	writeMultiBlockCustomSolution("out/xxixi.q",&gm,t.xxixi);
	writeMultiBlockCustomSolution("out/yxixi.q",&gm,t.yxixi);
	//writeMultiBlockGrid(&gm,"out/gridm.xyz");
	printf("Total new blocks: %d\n",gm.tb);
	printf("xd1: %d\n",gm.xdim[1][1]);
	printf("side: %d\n",m.sides[1]);
	int i,b;
	int**the=m.merges;
	for(i=0;i<m.tm;i++){
		printf("Merge %d: ",i);
		for(b=0;b<m.nb[i];b++){
			printf("%d ",the[i][b]);
		}
		printf("\n");
	}
	printf("Total merges: %d\n",m.tm);
	printf("Number of blocks in each merge: %d %d\n",m.nb[0],m.nb[1]);

	printf("Mask: %d %d %d %d %d %d %d %d\n",
			m.mask[0],m.mask[1],m.mask[2],m.mask[3],m.mask[4],m.mask[5],m.mask[6],m.mask[7]);

	printf("Total new blocks: %d\n",gm2.tb);
	printf("Key Test: %d %d %d %d\n",m.keynb[0],m.keynb[1],m.keynb[2],m.keynb[3]);
	printf("Key Test: %d %d %d\n",m.key[0][0],m.key[0][1],m.key[0][2]);
	printf("Key Test: %d %d %d\n",m.key[1][0],m.key[1][1],m.key[1][2]);
	printf("Sides: %d %d\n",m.sides[0],m.sides[1]);
	int i,j;
	for(j=0;j<10;j++){
		for(i=0;i<23;i++){
			printf("%g\t",gm2.x[1][i+j*gm2.xdim[1][0]]);
		}
		printf("\n");
	}
	int i,b;
	int**the=m.key;
	for(i=0;i<gm.tb;i++){
		printf("New Block %d: ",i);
		for(b=0;b<m.keynb[i];b++){
			printf("%d ",the[i][b]);
		}
		printf("\n");
	}
	*/

	printf("Finished!");
	return EXIT_SUCCESS;
}
