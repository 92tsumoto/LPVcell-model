#include <string.h>
#include "syspara.h"

void input_para(FILE *fpin)
{
	int i,ii;
	var.x0  = (double**)malloc2d(sizeof(double*),NUM,NN);

	fscanf(fpin,"%lf",&var.BCL);
	fscanf(fpin,"%lf",&var.ndis);
	fscanf(fpin,"%lf",&var.Istim_base);
	fscanf(fpin,"%lf",&Ligandtot);
	//fscanf(fpin,"%lf",&conIP3_inf);
	fscanf(fpin,"%lf",&conIP3);
	fscanf(fpin,"%lf",&scfIKb);
	fscanf(fpin,"%d",&var.datas);
	for (ii = 0; ii < var.datas; ii++){
		fscanf(fpin, "%d", &var.line_wid[ii]);
		for (i=0;i<NN;i++){
			fscanf(fpin,"%lf",&var.x0[ii][i]);
			//printf("%d,here\n",i);
			printf("x[%d]=%16.15e\n",i,var.x0[ii][i]);
		}
		fscanf(fpin, "%lf", &var.tsign[ii]);
		fscanf(fpin, "%lf", &var.tend[ii]);
		printf("tsign[%d]=%e\n",ii,var.tend[ii]);
	}
	fscanf(fpin,"%d",&var.l);
	fscanf(fpin,"%d",&var.m);
	fscanf(fpin,"%d",&var.pflag);
	fscanf(fpin,"%d",&var.write);
	fscanf(fpin,"%d",&var.write0);
	fscanf(fpin,"%d",&var.debug);

	xddp.line_wid = 1;
	xddp.line_att = SOLID;
	for (i=0;i<4;i++){
		fscanf(fpin,"%d",&xddp.sv[i]);
	}
	for (i=0;i<4;i++){
		fscanf(fpin,"%lf",&xddp.sc[i]);
	}
	for (i=0;i<2;i++){
		fscanf(fpin,"%lf",&xddp.la[i]);
	}
	for (i=0;i<2;i++){
		fscanf(fpin,"%d",&xddp.lav[i]);
	}
	for (i=0;i<2;i++){
		fscanf(fpin,"%d",&xddp.fi[i]);
	}
	fscanf(fpin,"%d",&xddp.kse);
	fscanf(fpin,"%s",xddp.font);
	for (i=0;i<4;i++){
		fscanf(fpin,"%d",&xddp.arrowtail[i]);
	}
	fscanf(fpin,"%d",&xddp.flag);
	printf("out input_para\n");
}

