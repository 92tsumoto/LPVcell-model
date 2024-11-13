/* produced by Tsumoto. K 2008.10.27 */
/* produced by Tsumoto. K 2022.11.28 */
/* produced by Tsumoto. K 2024.09.29 */

#include <string.h>
#include <stdlib.h>
#include "syspara.h"

FILE *fpin; 
FILE *fp0,*fp1,*fp2,*fp3,*fp4,*fp5;
FILE *fp6,*fp7,*fp8,*fp9,*fp10;
int mode = 1;
int P = 2;
int beats = 1;

typedef double Number;

int main(int argc, char **argv)
{
	int i,w;
	int ii=0;
	double x[NN];
	double t = 0.0;
	double tt;
	double time=0.0;
	double ttime=0.0;
	double h;
	double v_old,dvdt_new;
	char *tmpname;
	char cmd[BUFSIZ];
	double tend;

/* Action Potential Duration and Max. Info */
	double *vmax ; // Max. Voltage (mV)
	double *dvdtmax ; // Max. dv/dt (mV/ms)
	double *apd; // Action Potential Duration
	double *toneapd; // Time of dv/dt Max.
	double *ttwoapd; // Time of 90% Repolarization
	double *rmbp; // Resting Membrane Potential
	double *nair; // Intracellular Na At Rest
	double *cair; // Intracellular Ca At Rest
	double *kir ; // Intracellular K At Rest

	vmax=(Number *)calloc(beats,sizeof(Number));
	dvdtmax=(Number *)calloc(beats,sizeof(Number));
	apd=(Number *)calloc(beats,sizeof(Number));
	toneapd=(Number *)calloc(beats,sizeof(Number));
	ttwoapd=(Number *)calloc(beats,sizeof(Number));
	rmbp=(Number *)calloc(beats,sizeof(Number));
	nair=(Number *)calloc(beats,sizeof(Number));
	cair=(Number *)calloc(beats,sizeof(Number));
	kir=(Number *)calloc(beats,sizeof(Number));
	if(vmax==NULL || dvdtmax==NULL || apd==NULL || toneapd==NULL || ttwoapd==NULL 
		|| rmbp==NULL || nair==NULL || cair==NULL || kir==NULL
		) exit(1);

	tmpname = "temp";

	sprintf(cmd, "/usr/bin/cpp -P %s > %s", argv[1],tmpname);
	if(system(cmd) == -1){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if((fpin=fopen(tmpname,"r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if ((fp1 = fopen("para.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp2 = fopen("data.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp3 = fopen("nstate.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}

// parameter inputs
	input_para(fpin);

	if (var.write){
		if ((fp0 = fopen(argv[2],"w"))==NULL){
			fprintf(stderr, "%s cannot open.\n",argv[2]);
			exit(-1);
		}
	}
	if (var.write0){
		if ((fp4 = fopen("ikr.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp5 = fopen("ina.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp6 = fopen("ical.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp7 = fopen("incx.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp8 = fopen("inak.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp9 = fopen("ikur.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
		if ((fp10 = fopen("ikto.out","w")) == NULL){
			printf("Can't open File\n");
			exit(1);
		}
	//	if ((fp9 = fopen("jrel.out","w")) == NULL){
	//		printf("Can't open File\n");
	//		exit(1);
	//	}
	}

	xhplot(WINDOW, 700.0, 700.0, WHITE);
	xhplot(DIRECT, 0.0, 0.0, WHITE);

	for (ii = 0; ii < var.datas; ii++){
		long j;
		time = 0.0;
		tend = var.tend[ii];
		for (i = 0; i < NN; i++){ 
			x[i] = var.x0[ii][i];
		}

		//tt = var.ndis*(double)var.m;
		h = 1.0 / (double)var.m;
		//h = 2.0*M_PI / (double)var.m;
		h *= var.tsign[ii];

		xddp.line_wid = var.line_wid[ii];
		xhplot(LINEATT,0,0,WHITE);

		// initial values input.
		val_consts(fp1);
	
		// initial values input.
		initial_mem();
		printf("exit memory initialization\n");

		printf("Istim=%lf\n",var.Istim_base);

		// Tablize exp functions.	
		printf("start tablization\n");
		make_ExpTable();
		printf("finished tablization\n");

		// Initialization time
		var.beat = 0;

		//tt = var.ndis*(double)var.m*var.BCL;
		tt = (double)var.m*var.BCL;
		//ttt = (1.0-var.ndis)*(double)var.m*var.BCL;
		printf("tt=%lf,time=%lf\n",tt,time);
		var.Istim = var.Istim_base;

		while (1){
			eventloop(fp1,&mode,&P,x);
			ttime=var.beat*var.BCL;

			for (j = 0; j< (int)tt; j++){
				t = h*(double)j;
				v_old = x[0];
			/*if(var.beat == 0){
				if(time-(ttime+20.0) >= 0.0 && time-(ttime+20.0) < 3.0 ){
				//if(time-(ttime+20.0) >= 0.0 && time-(ttime+20.0) < 180.0 ){
					var.Istim = var.Istim_base;
				} else {
					var.Istim = 0.0;
				}
			} else { var.Istim = 0.0; }*/
				runge(NN,h,x,t);
				//eular(NN,h,x,t);
				dvdt_new = (x[0]-v_old)/h;
				if (var.pflag) orbit(&mode,x,dvdt_new);
				if (var.pswitch==1){
				//if (1){
					//printf("data out\n");
					data_out(fp2,time,x);
					if(var.write0){
						current(fp4,fp5,fp6,fp7,fp8,fp9,fp10,time,x);
					}
				}
				time += h;
				//ttime=time;
			}

			fprintf(fp3,"#beats=%d\n",var.beat);
			for(w=0;w<NN;w++){
				fprintf(fp3,"%16.15e\n",x[w]);
			}
			fflush(fp3);

			printf("%d %lf ",var.beat,time);
			for(w=0;w<NN;w++){
				if(w!=NN-1){
					printf("%10.9lf ",x[w]);
				} else {
					printf("%10.9lf %lf %lf\n",x[w],Ligand,conIP3_inf);
				}
			}
			draw_p(&mode,P,x,dvdt_new);
			mouse(&mode,x,dvdt_new);
			if (fabs(time) > tend &&  tend != 0.0) break;
			var.beat++;

		} // end for while loop

	} // end for ii-loop

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	if(var.write0){
		fclose(fp4);fclose(fp5);fclose(fp6);fclose(fp7);fclose(fp8);
	}
	free(vmax);free(dvdtmax);free(apd);free(toneapd);free(ttwoapd);
	free(rmbp);free(nair);free(cair);free(kir);

}

