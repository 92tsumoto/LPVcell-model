#include "syspara.h"

void make_ExpTable()
{

	//int vindex,kiindex;
	int vindex;
	//double v,ki;
	double v;
    
	for(vindex=0;vindex<VNMAX;vindex++){

        v = (double)vindex/dvm-Emax;
		
		/** for ina **/
		TkC2O[vindex] = 0.5/(0.0025*exp(-(v-10.0)/8.0)+0.15*exp(-(v-10.0)/100.0));
		TkOC[vindex] = 0.5/(30.0*exp((v-10.0)/12.0)+0.53*exp((v-10.0)/50.0));
		TkOI2[vindex] = 0.1/(0.0433*exp(-(v-10.0)/27.0) + 0.34*exp(-(v-10.0)/2000.0));	 // Umehara et al., PVC model

		TkIsb[vindex] = 1.0/(150000.0*exp((v-15.0)/10.0)+25000.0*exp((v-15.0)/16.0)); 	// Umehara et al., PVC model
		TkIsf[vindex] = 1.0/(0.008*exp(-(v-15.0)/9.9)+4.0*exp(-(v-15.0)/45.0));		// Umehara et al., PVC model
		TfC[vindex] = 1.0/(1.0+exp(-(v+38.0)/7.0));

		// iKur
		Tass[vindex] = 1.0/(1.0+exp(-(v+22.5)/7.7));
		Tiss[vindex] = 1.0/(1.0+exp((v+45.2)/5.7));
		Ttau_aur[vindex] = 0.493*exp(-0.0629*v) + 2.058;
		Ttau_iur[vindex] = 1200.0-170.0/(1.0+exp((v+45.2)/5.7));

		// ito
		// Ca independent transient outward K current,r-gate
		Tpr_inf[vindex] = 1.0/(1.0+exp((v+10.6) / -11.42));
		Ttau_pr[vindex] = 1000.0/(45.16*exp(0.03577*(v+50.0)) + 98.9*exp(-0.1*(v+38.0)));
		// Ca independent transient outward K current,s-gate: inactivation gate
		Ts_inf[vindex] = 1.0/(1.0+exp((v+45.3)/6.8841));
		Ttau_s[vindex] = 1000*(0.35*exp(-(v+70.0)*(v+70)/15.0/15.0)+0.035);
		// Ca independent transient outward K current, s-slow-gate: inactivation gate
		Tss_inf[vindex] = 1.0/(1.0+exp((v+45.3)/6.8841));
		Ttau_ss[vindex] = 1000.0*10.0*(3.7*exp(-(v+70.0)*(v+70.0)/30.0/30.0)+0.035); // decelerated by 10000 'use original model

		// for ical
		TfVmACT[vindex] = 1.0/(3.734*exp(-(v + VshiftLCCact)/8.5) + 0.35*exp(-(v + VshiftLCCact)/3500.0));
		TexpF[vindex] = exp(-v/RTon2F);
		Tkooco_iz[vindex] = 1.0 / (4.65*exp((v + VshiftLCCinact)/15.0) + 1.363*exp((v + VshiftLCCinact)/100.0));
		Tkccco_iz[vindex] = 1.0/(8084.0*exp(v/10.0) + 158.0*exp(v/1000.0)) + 1.0/(134736.0*exp(-v/5.0) + 337.0*exp(-v/2000.0));
		Tkooco_blk[vindex] = 1.0 / (4.65*exp((v + VshiftLCCinact)/15.0) + 1.363*exp((v + VshiftLCCinact)/100.0));
		Tkccco_blk[vindex] = 1.0/(8084.0*exp(v/10.0) + 158.0*exp(v/1000.0)) + 1.0/(134736.0*exp(-v/5.0) + 337.0*exp(-v/2000.0));

        //'**********  IKr modified by U. for PVC at the initial stage of study *********************
        // the fast component of activation
        Tay1[vindex] = 1.0 / (35.0 * exp(-v / 10.5) + 75.0 * exp(-v / 100.0));
        Tby1[vindex] = 1.0 / (470.0 * exp(v / 8.3) + 220.0 * exp(v / 29.0));
        // the slow component of activation
        Tay2[vindex] = 1.0 / (350.0 * exp(-v / 10.5) + 300.0 * exp(-v / 100.0));
        Tby2[vindex] = 1.0 / (1850.0 * exp(v / 8.3) + 2200.0 * exp(v / 29.0));
        // inactivation gate (fast gate)
        Tay3[vindex] = 1.0 / (0.015 * exp(v / 6.0) + 7.0 * exp(v / 60.0));
        Tby3[vindex] = 1.0 / (0.114 * exp(-v / 9.2) + 2.3 * exp(-v / 1000.0));

        
		// iClh
        Ta1_Clh[vindex] = 0.001 * 0.064 * exp(-0.028 * v);
        Ta2_Clh[vindex] = 0.001 * 0.8 * exp(-0.02 * v);
        Tb1_Clh[vindex] = 0.001 * 7.5 * exp(0.018 * v);
        Tb2_Clh[vindex] = 0.001 * 1.0 / (0.0006 * exp(-v / 15.6) + 0.0594 * exp(-v / 126));
        
		// iNCX
		Tk1_NCX[vindex] = exp(0.32*v/RTonF);	// exp(partition * x[0]/RTonF)
		Tk2_NCX[vindex] = exp(-0.68*v/RTonF);	// exp((partition - 1.0)*x[0]/RTonF)

		// RyR
		Tkvoc[vindex] = 1.0/(4.65*exp((v + VshiftLCCinact)/15.0) + 1.363*exp((v + VshiftLCCinact)/100.0));
		Tkcaco[vindex] = 1.0/(8084.0*exp(v/10.0)+158.0*exp(v/1000.0))+1.0/(134736.0*exp(-v/5.0) + 337.0*exp(-v/2000.0));
	}

}
