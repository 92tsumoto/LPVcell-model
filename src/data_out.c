#include "syspara.h"

void data_out(FILE *fp2, double t, double u[])
{
	int i;

	printf("%lf ",t);
	fprintf(fp2,"%lf ",t);
	for(i=0;i<NN;i++){
		printf("%lf ",u[i]);
		fprintf(fp2,"%lf ",u[i]);
	}
	//fprintf(fp2,"%lf %lf\n",t,u[0]);
	fprintf(fp2,"\n");
	printf("\n");
	//printf("output data end\n");

}

void current(FILE *fp4, FILE *fp5, FILE *fp6, FILE *fp7, FILE *fp8, FILE *fp9, FILE *fp10, double t, double u[])
{

	current_ikr(fp4,t,u);
	current_ina(fp5,t,u);
	current_ical(fp6,t,u);
	current_inaca(fp7,t,u);
	current_inak(fp8,t,u);
	current_ikur(fp9,t,u);
	current_ikto(fp10,t,u);
	//printf("t=%lf\n",t);

}

// L-type calcium current
void current_ical(FILE *fp6, double time, double p[])
{

	fprintf(fp6,"%lf %lf %lf\n",time,ICaLCa_iz,ICaLCa_blk);

}

// Rapidly Activating Potassium Current 
void current_ikr (FILE *fp4, double time, double p[])
{
	fprintf(fp4,"%lf %lf\n",time,IKr);

}

void current_ikur (FILE *fp9, double time, double p[])
{
	fprintf(fp9,"%lf %lf\n",time,IKur);

}

void current_ikto (FILE *fp10, double time, double p[])
{
	fprintf(fp10,"%lf %lf\n",time,IKto);

}

// Slowly Activating Potassium Current 
void current_ina (FILE *fp5, double time, double p[])
{
	
	fprintf(fp5,"%lf %lf %lf %lf\n",time,INaT_Na,INaT_K,INaT);

}

// Sodium-Calcium Exchanger V-S

void current_inaca (FILE *fp7, double time, double p[])
{
	fprintf(fp7,"%lf %lf %lf %lf\n",time,totINCX,INCXNa_jnc,INCXCa_jnc);

}

// Sodium-Potassium Pump

void current_inak (FILE *fp8, double time, double p[])
{
	fprintf(fp8,"%lf %lf %lf %lf\n",time,INaK,INaK_Na,INaK_K);

}
