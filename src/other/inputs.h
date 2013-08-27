void numericsDict(Numerics* n,char *var, double* val){
	int matchlimit=4;
	if(strncasecmp(var, "useControlTerms", matchlimit)==0){ n->useControlTerms=val[0]; }
	if(strncasecmp(var, "writeResiduals", matchlimit)==0){ n->writeResiduals=val[0]; }
	if(strncasecmp(var, "tol", matchlimit)==0){ n->tol=val[0]; }
	if(strncasecmp(var, "difftol", matchlimit)==0){ n->difftol=val[0]; }
	if(strncasecmp(var, "cornertol", matchlimit)==0){ n->cornertol=val[0]; }
	if(strncasecmp(var, "spikelim", matchlimit)==0){ n->spikelim=val[0]; }
}
void readNumerics(Numerics* n,char* fn){
	int maxchars=200;
	char line[maxchars];
	char var[maxchars];
	double val[2];
	FILE * file = fopen (fn, "rt");
	int success;
	//Read total number of blocks.
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf(line, "%s %lf %lf",var,&val[0],&val[1]);
		if(success>0){ numericsDict(n,var,val); }
	}
	fclose(file);
}
