char **charmalloc(int strings,int length){
	char **chararray=(char **)malloc(sizeof(char *)*strings);
	int i;
	for(i=0;i<strings;i++){
		chararray[i]=malloc(sizeof(char)*length);
	}
	return chararray;
}
int **mbinfo(int blocks){
	int b;
	int **mbi=(int **)malloc(sizeof(int *)*blocks);
	for(b=0;b<=blocks;b++){
		mbi[b]=(int *)malloc(sizeof(int)*4);
	}
	return mbi;
}
void mbcopy(Grid*g,double**from,double**to){
	int b,xi,xd3;
	for(b=0;b<g->tb;b++){
		xd3=g->xdim[b][3];
		for(xi=0;xi<xd3;xi++){
			to[b][xi]=from[b][xi];
		}
	}
}
void initv(double * vec, int len, double val){
	int i;
	for (i=0;i<len;i++){
		vec[i]=val;
	}
}
