#include "syspara.h"

void comp_potential(double x[])
{
	GHK_Na = zna*x[0]/RTonF*(x[na]i-x[na]o*exp(-zna*x[0]/RTonF))/(1.0 - exp(-zna*x[0]/RTonF));
	GHK_K  = zk *x[0]/RTonF*(x[k]i -x[k]o *exp(-zk *x[0]/RTonF))/(1.0 - exp(-zna*x[0]/RTonF));
	GHK_Ca = zca*x[0]/RTonF*(x[ca]i-x[ca]o*exp(-zca*x[0]/RTon2F))/(1.0 - exp(-zna*x[0]/RTonF));
}

void comp_ina(double x[])
{
	double vm;
	int iV = 0;
	double V1,V2,d1,d2;
	vm = x[0];
	pC_NaT = x[1];
	PO_NaT = x[2];
	pI2_NaT = x[3];
	pIs_NaT = x[4];

	V1 = (vm + Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

// Transient component (fast?)
	fc_Na = Tfc_Na[iV]*d2 + Tfc[iV+1]*d1;
	k_C2O = Tk_C2O[iV]*d2 + Tk_C2O[iV+1]*d1;
	k_OC  = Tk_OC[iV]*d2  + Tk_OC[iV+1]*d1;
	k_OI2 = Tk_OI2[iV]*d2 + Tk_OI2[iV+1]*d1;
	k_Isb = Tk_Isb[iV]*d2 + Tk_Isb[iV+1]*d1;
	k_Isf = Tk_Isf[iV]*d2 + Tk_Isf[iV+1]*d1;
	
	k_C2I2 = 0.5/(1.0 + k_I2O*k_OC/k_OI2/k_C2O);
	k_I2C = 0.5 - k_C2I2;

	ina.fast = (1.0 - fL)*Pna*(GHK_Na + 0.1*GHK_K)*x[1];

// Late component
	k_OI1 = k_OI2;
	k_I1C = k_I2C;
	k_C2I1 = k_C2I2;

	ina.late = fL*PNa*(GHK_Na+0.1*GHK_K)*x[6];
	
	ina.total = ina.fast + ina.late;
}

