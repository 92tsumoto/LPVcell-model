#include "syspara.h"

void comp_potential(double x[])
{
	double vm;	
	double Cafree_jnc,Cafree_iz,Cafree_blk;
	double Nai,Ki;

	vm         = x[0];
	Cafree_jnc = x[39];
	Cafree_iz  = x[40];
	Cafree_blk = x[41];
	Nai        = x[44];
	Ki         = x[45];

	Ena     = RTonF*log(Nao/Nai);		
	Ek      = RTonF*log(Ko/Ki);			
	ECa_jnc = RTon2F*log(Cao/Cafree_jnc);	
	ECa_iz  = RTon2F*log(Cao/Cafree_iz );	
	ECa_blk = RTon2F*log(Cao/Cafree_blk);	
	//ECa_LR  = RTon2F*log(Cao/Ca_ndLR);
	//ECa_L0  = RTon2F*log(Cao/Ca_ndL0);

	GHK_Na = zna*vm*F/R/T*(Nai - Nao*exp(-zna*vm*F/R/T))/(1.0 - exp(-zna*vm*F/R/T));
	GHK_K  = zk *vm*F/R/T*(Ki  - Ko *exp(-zk *vm*F/R/T))/(1.0 - exp(-zk *vm*F/R/T));
	GHK_Ca_jnc = zca*vm/RTonF*(Cafree_jnc  - Cao*exp(-zca*vm/RTonF))/(1.0 - exp(-zca*vm/RTonF));
	GHK_Ca_iz  = zca*vm/RTonF*(Cafree_iz   - Cao*exp(-zca*vm/RTonF))/(1.0 - exp(-zca*vm/RTonF));
	GHK_Ca_blk = zca*vm/RTonF*(Cafree_blk  - Cao*exp(-zca*vm/RTonF))/(1.0 - exp(-zca*vm/RTonF));
	//GHK_Ca_LR  = zca*vm/RTonF*(Ca_ndLR - Cao*exp(-zca*vm/RTonF))/(1.0 - exp(-zca*vm/RTonF));
	//GHK_Ca_L0  = zca*vm/RTonF*(Ca_ndL0 - Cao*exp(-zca*vm/RTonF))/(1.0 - exp(-zca*vm/RTonF));
}

void comp_ina (double x[])
{

	double vm,O_TM,I2_TM,Is_TM;

	vm     = x[0];
    O_TM   = x[1];
    I2_TM  = x[2];
    Is_TM  = x[3];

	kC2O  = 0.5/(0.0025*exp(-(vm-10.0)/8.0)+0.15*exp(-(vm-10.0)/100.0));
	kOC   = 0.5/(30.0*exp((vm-10.0)/12.0)+0.53*exp((vm-10.0)/50.0));
	kOI2  = 0.1/(0.0433*exp(-(vm-10.0)/27.0) + 0.34*exp(-(vm-10.0)/2000.0));   // Umehara et al., PVC model
	//kI2O  = 0.0001312;

	kC2I2 = 0.5/(1.0 + kI2O*kOC/kOI2/kC2O);
	kI2C  = 0.5 - kC2I2;
	
	kIsb  = 1.0/(150000.0*exp((vm-15.0)/10.0)+25000.0*exp((vm-15.0)/16.0));
	kIsf  = 1.0/(0.008*exp(-(vm-15.0)/9.9)+4.0*exp(-(vm-15.0)/45.0));

	fC  = 1.0/(1.0+exp(-(vm+38.0)/7.0));

	C_TM = 1.0 - Is_TM - O_TM - I2_TM;

	INaT_Na = rate_INa *PNa*GHK_Na*O_TM;
	INaT_K  = rate_INa *PNa*rPK_Ina*GHK_K*O_TM;
	INaT = INaT_Na + INaT_K;

}

// IKur Vladimir E. Bondarenko, et al. mouse ventricular myocytes model
void comp_ikur (double x[])
{

	double vm,a_Kur,i_Kur;
	
	vm     = x[0];
    a_Kur  = x[4];
    i_Kur  = x[5];

	ass      = 1.0/(1.0+exp(-(vm+22.5)/7.7));
	iss      = 1.0/(1.0+exp((vm+45.2)/5.7));
	tau_aur  = 0.493*exp(-0.0629*vm) + 2.058;
	tau_iur  = 1200.0-170.0/(1.0+exp((vm+45.2)/5.7));

 	IKur = rate_Kur*GKur*a_Kur*i_Kur*(vm - Ek);
}

// Ito Pandit et al.
void comp_ito (double x[])
{

	double vm,pr,s_slow;
	
	vm		= x[0];
    pr		= x[6];
    s_slow	= x[7];

    // Ca independent transient outward K current,r-gate
	pr_inf	= 1.0/(1.0+exp( (vm + 10.6) / -11.42) );
	tau_pr	= 1000.0/(45.16*exp(0.03577*(vm + 50.0)) + 98.9*exp(-0.1*(vm+38.0)));
    // Ca independent transient outward K current,s-gate: inactivation gate
    //s_inf   = 1.0/(1.0+exp((vm + 45.3)/6.8841));
    //tau_s   = 1000.0*(0.35*exp(-(vm + 70.0)*(vm + 70)/15.0/15.0) + 0.035);
    // Ca independent transient outward K current, s-slow-gate: inactivation gate
    ss_inf  = 1.0/(1.0+exp((vm + 45.3)/6.8841));
    tau_ss  = 1000.0*10.0*(3.7*exp(-(vm + 70.0)*(vm + 70.0)/30.0/30.0)+0.035); // decelerated by 10000 'use original model

	IKto = rate_Ito*GKto*pr*s_slow*(vm - Ek);

}


// L-type calcium current
void comp_ical(double x[])
//void comp_ical_iz(double x[])
{
	double Ca_loc_iz,Ca_loc_blk; // local [Ca] at the Ca binding site (calmodulin tethered to LCC)
	double Cafree_iz,Cafree_blk;
	double Yco_iz,Yoo_iz,Yoc_iz;
	double Yco_blk,Yoo_blk,Yoc_blk;
	double vm,fVmAct,expF;
	double product_iz,product_blk;
	double perNa_ICaL,perK_ICaL;
	
	vm         =  x[0];
	Yco_iz     =  x[8]; // Yco_iz
	Yoc_iz     =  x[9]; // Yoc_iz
	Yoo_iz     = x[10]; // Yoo_iz
	Yco_blk    = x[11];
	Yoc_blk    = x[12];
	Yoo_blk    = x[13];
	Cafree_blk = x[41];
	Cafree_iz  = x[40];

	expF	= exp(-vm/RTon2F);
	//fVmAct  = 1.0/(3.734*exp(-(vm + VshiftLCCact)/8.5) + 0.35*exp(-(vm + VshiftLCCact)/3500.0));
	fVmAct  = 1.0/(3.734*exp(-vm/8.5) + 0.35*exp(-vm/3500.0));
	Ca_loc_iz  = (Cafree_iz  + JLCC/gD_free*Cao*(vm/RTon2F*expF)/(1.0 - expF))/(1.0 + JLCC/gD_free*vm/RTon2F/(1.0 - expF));
	Ca_loc_blk = (Cafree_blk + JLCC/gD_free*Cao*(vm/RTon2F*expF)/(1.0 - expF))/(1.0 + JLCC/gD_free*vm/RTon2F/(1.0 - expF));

/*	if(fabs(vm) < 0.00001){
		//Ca_loc = -(Cafree_blk - JLCC/gD_free*Cao/RTon2F)/(JLCC/gD_free/RTon2F - 1.0);
		Ca_loc = (Cafree_iz + JLCC*Cao/gD_free)/(1 + JLCC/gD_free);
	} else {
		Ca_loc = (Cafree_iz + JLCC/gD_free*Cao*(vm/RTon2F*expF)/(1.0 - expF))/(1.0 + JLCC/gD_free*vm/RTon2F/(1.0 - expF));
	}
*/
	kcooo_iz= fVmAct;
	kooco_iz= 1.0 / (4.65*exp((vm + VshiftLCCinact)/15.0) + 1.363*exp((vm + VshiftLCCinact)/100.0));
	kccoc_iz= fVmAct;
	koccc_iz= 1.0 / (4.65*exp((vm + VshiftLCCinact)/15.0) + 1.363*exp((vm + VshiftLCCinact)/100.0));

	kccco_iz= 1.0/(8084.0*exp(vm/10.0) + 158.0*exp(vm/1000.0)) + 1.0/(134736.0*exp(-vm/5.0) + 337.0*exp(-vm/2000.0));
	kcocc_iz = (Cafree_iz/KL)*fVmAct/TL;
	kocoo_iz= 1.0/(8084.0*exp(vm/10.0) + 158.0*exp(vm/1000.0)) + 1.0/(134736.0*exp(-vm/5.0) + 337.0*exp(-vm/2000.0));
	koooc_iz = (Ca_loc_iz/KL)*fVmAct/TL;

	kcooo_blk= fVmAct;
	kooco_blk= 1.0 / (4.65*exp((vm + VshiftLCCinact)/15.0) + 1.363*exp((vm + VshiftLCCinact)/100.0));
	kccoc_blk= fVmAct;
	koccc_blk= 1.0 / (4.65*exp((vm + VshiftLCCinact)/15.0) + 1.363*exp((vm + VshiftLCCinact)/100.0));

	kccco_blk= 1.0/(8084.0*exp(vm/10.0) + 158.0*exp(vm/1000.0)) + 1.0/(134736.0*exp(-vm/5.0) + 337.0*exp(-vm/2000.0));
	kcocc_blk = (Cafree_blk/KL)*fVmAct/TL;
	kocoo_blk= 1.0/(8084.0*exp(vm/10.0) + 158.0*exp(vm/1000.0)) + 1.0/(134736.0*exp(-vm/5.0) + 337.0*exp(-vm/2000.0));
	koooc_blk = (Ca_loc_blk/KL)*fVmAct/TL;
	
	Ycc_iz  = 1.0 - (Yco_iz  + Yoc_iz  + Yoo_iz);
	Ycc_blk = 1.0 - (Yco_blk + Yoc_blk + Yoo_blk);

//printf("kccco_iz=%lf,kooco_iz=%lf,kcocc_iz=%lf,kcooo_iz=%lf,Cafree_iz=%lf\n",kccco_iz,kooco_iz,kcocc_iz,kcooo_iz,x[42]);
//  afLCC = baseCaL + deltaCaL*y_PKA;
//	afLCC = 1.0 + 1.3*y_PKA;
	
	product_iz  = Yoo_iz  * ATPfactor;
	product_blk = Yoo_blk * ATPfactor;
	tPCaL = (oPCaL * scf_LCC) * afLCC;

	perNa_ICaL = 0.0000185*tPCaL;
	perK_ICaL  = 0.000367*tPCaL;

	ICaLCa_iz = Frc_CaL_iz * tPCaL * GHK_Ca_iz * product_iz;	// (pA/pF): amplitude of Ca component
	ICaLNa_iz = Frc_CaL_iz * perNa_ICaL * GHK_Na * product_iz;	// (pA/pF), amplitude of  Na component
	ICaLK_iz  = Frc_CaL_iz * perK_ICaL * GHK_K * product_iz;	// (pA/pF), amplitude of K component

	ICaLCa_blk = Frc_CaL_blk * tPCaL * GHK_Ca_blk * product_blk;	// (pA/pF): amplitude of Ca component
	ICaLNa_blk = Frc_CaL_blk * perNa_ICaL * GHK_Na * product_blk;	// (pA/pF), amplitude of  Na component
	ICaLK_blk  = Frc_CaL_blk * perK_ICaL * GHK_K * product_blk;	// (pA/pF), amplitude of K component

}

// Rapidly Activating Potassium Current 
void comp_ikr (double x[])
{

	double vm,y1,y2,y3,pO;

	vm = x[0];
	y1 = x[14];
	y2 = x[15];
	y3 = x[16];

	ay1	= 1.0 / (35.0 * exp(-vm / 10.5) + 75.0 * exp(-vm / 100.0));
	by1	= 1.0 / (470.0 * exp(vm / 8.3) + 220.0 * exp( vm / 29.0));
    // the slow component of activation
	ay2	= 1.0 / (350.0 * exp(-vm / 10.5) + 300.0 * exp(-vm / 100.0));
	by2	= 1.0 / (1850.0 * exp(vm / 8.3) + 2200.0 * exp( vm / 29.0));
    // inactivation gate (fast gate)
	ay3	= 1.0 / (0.015 * exp( vm / 6.0) + 7.0 * exp( vm / 60.0));
	by3	= 1.0 / (0.114 * exp(-vm / 9.2) + 2.3 * exp(-vm / 1000.0));

	pO = (0.6*y1 + 0.4*y2)*y3;

	IKr = rate_Kr*gKr*qKr*pO*(vm - Ek);

}

// Inward rectifier potassium current (Ik1)
void comp_ik1 (double x[])
{
	//double vm,Pbspm,ssPbspm;
	double vm,Pbspm;
	double aMg,bMg;
	double fracO,fracB;
	double poMg1,poMg2;
	double pO_IK1,pO_mode1,pO_mode2;
	
	vm    = x[0];
	Pbspm = x[17];

	//'************** Mg-block *************************
	aMg = 12.0*exp(-0.025*(vm - Ek)); 
	bMg = 28.0*Mg_cyt*exp(0.025*(vm - Ek));
	
	fracO = aMg / (aMg + bMg);
	fracB = bMg / (aMg + bMg);

	poMg2 = 3.0 * fracO * fracB * fracB;
	poMg1 = 3.0 * fracO * fracO * fracB;
	po_Mggate = fracO * fracO * fracO;
        
	//'************** spermin-block *************************
	aSPM = 0.17 * exp(-0.07 * ((vm - Ek) + 8.0*Mg_cyt))/(1.0+0.01*exp(0.12*((vm - Ek) + 8.0*Mg_cyt)));
	bSPM = 0.28*SPM*exp(0.15*((vm - Ek) + 8.0*Mg_cyt))/(1.0+0.01*exp(0.13*((vm - Ek) + 8.0*Mg_cyt)));

	//f[18]:dPbspmdt; = bSPM*po_Mggate*(1.0-x[18]) - aSPM*x[18];     'time dependent gating
	//ssPbspm = (bSPM * po_Mggate) / (aSPM + bSPM*po_Mggate);	//  instantaneous equilibrium

	//'**************** IK1 amplitude **************************************************
	pO_mode1 = frac_mode_K1*(1.0 - Pbspm)*(po_Mggate + 2.0*poMg1/3.0 + 1.0*poMg2/3.0);     //'time-dependent gating
	//printf("m1=%9.8lf,pom0=%9.8lf,pom1=%9.8lf,pom2=%9.8lf\n",frac_mode_K1*(1.0 - Pbspm),po_Mggate,2.0*poMg1/3.0,1.0*poMg2/3.0);
	pO_mode2 = (1.0 - frac_mode_K1) / (1.0 + SPM / (40.0*exp(-(vm - Ek)/9.1)));

	pO_IK1 = (pO_mode1 + pO_mode2);

	//fK1 = qK1/(1.0 + exp(-(Ko - 2.2)/0.6));

	IK1 = GK1*fK1*pO_IK1*(vm-Ek);
	//printf("cond=%9.8lf,fk1=%9.8lf,pO=%9.8lf,driv=%9.8lf\n",GK1,fK1,pO_IK1,vm-Ek);
}

// OKamoto et al., JMCC, 2012,2014
// IClh: a specific cloride current 
void comp_iclh (double x[])
{

	double vm,c1,c2,O,Ecl;
	// double o_clh--> none;	//state variables {init: 0.0219499323546152}
	// double c1_clh-->x[19];	//state variables {init: 0.118714932500685}
	// double c2_clh--> x[20];	//state variables {init: 0.859335135183572}

	vm	= x[0];
	c1	= x[18];
	c2	= x[19];
	Ecl = -20.0; // just an assumption; ECl = -20 mV	

	a1_Clh = 0.001 * 0.064 * exp(-0.028 * vm);
    a2_Clh = 0.001 * 0.8 * exp(-0.02 * vm);
    b1_Clh = 0.001 * 7.5 * exp(0.018 * vm);
    b2_Clh = 0.001 * 1.0 / (0.0006 * exp(-vm / 15.6) + 0.0594 * exp(-vm / 126.0));

	O = 1.0 - c1 - c2;
        
	IClh = rate_Clh * GClh * O * (vm - Ecl);            // 'just an assumption; ECl = -20 mV

}


// Sodium-Calcium Exchanger V-S
void comp_INCX (double x[])
{
	
	double E1Na_jnc,E1Ca_jnc;
	double E1Na_iz,E1Ca_iz;
	double E1Na_blk,E1Ca_blk;
	double fCaina_jnc,fCaina_iz,fCaina_blk;
	double k1,k2;
	double vm;
	double I1NCX_jnc,I2NCX_jnc,E1NCX_jnc;
	double I1NCX_iz,I2NCX_iz,E1NCX_iz;
	double I1NCX_blk,I2NCX_blk,E1NCX_blk;
	double Nai,Cafree_jnc,Cafree_iz,Cafree_blk;
	double partition;

	vm	       = x[0];
	E1NCX_jnc  = x[18];
	I1NCX_jnc  = x[19];
	I2NCX_jnc  = x[20];
	E1NCX_iz   = x[21];
	I1NCX_iz   = x[22];
	I2NCX_iz   = x[23];
	E1NCX_blk  = x[24];
	I1NCX_blk  = x[25];
	I2NCX_blk  = x[26];
	Cafree_jnc = x[39];
	Cafree_iz  = x[40];
	Cafree_blk = x[41];
	Nai        = x[44];

	partition = 0.32;
	k1 = exp(      partition*vm/RTonF);
	k2 = exp((partition-1.0)*vm/RTonF);

	// jnc NCX
	//E2Na = 1.0 / (1.0 + pow((KmNao / Nao),3.0) * (1.0 + Cao / KmCao));	// Na binding probability
	//E2Ca = 1.0 / (1.0 + (KmCao / Cao) * (1.0 + pow((Nao / KmNao),3.0)));	// Ca binding probability
	E1Na_jnc = 1.0 / (1.0 + pow((KmNai/Nai),3.0) * (1.0 + Cafree_jnc / KmCai));	// Na binding probability
	E1Ca_jnc = 1.0 / (1.0 + (KmCai/Cafree_jnc) * (1.0 + pow((Nai/KmNai),3.0)));	// Ca binding probability
	fCaina_jnc = Cafree_jnc/(Cafree_jnc + Kmact);	// Ca-dependency of the inactivation step

	// I1 inactivation gate (Na-dependent)
	alpha1_jnc = (E1Na_jnc*(fCaina_jnc * a1on + (1.0 - fCaina_jnc)*a1off));// *fixzero_slowVariable;	// rate from E1 to I1
	beta1_jnc  = (fCaina_jnc*b1on + (1.0 - fCaina_jnc)*b1off); // *fixzero_slowVariable;	// rate from I1 to E1
	// I2 inactivation gate (Ca-dependent)
	alpha2_jnc = (fCaina_jnc*a2on + (1.0 - fCaina_jnc)*a2off);	// rate from E1 to I2
	beta2_jnc  = (fCaina_jnc*b2on + (1.0 - fCaina_jnc)*b2off);	// rate from I2 to E1
	alphaE_jnc = k2*E2Na + k4*E2Ca;	// overall transition rate from E2 to E1
	betaE_jnc  = k1*E1Na_jnc + k3*E1Ca_jnc;	// overall transition rate from E1 to E2

	E2NCX_jnc = 1.0 - E1NCX_jnc - I1NCX_jnc - I2NCX_jnc;
	INCX_jnc = Frc_NCX_jnc*maxINCX*( k1*E1Na_jnc*E1NCX_jnc - k2*E2Na*E2NCX_jnc );
	INCXNa_jnc = 3.0 * INCX_jnc;
	INCXCa_jnc = -2.0 * INCX_jnc;
	
	// iz NCX
	E1Na_iz = 1.0 / (1.0 + pow((KmNai/Nai),3.0) * (1.0 + Cafree_iz / KmCai));	// Na binding probability
	E1Ca_iz = 1.0 / (1.0 + (KmCai/Cafree_iz) * (1.0 + pow((Nai/KmNai),3.0)));	// Ca binding probability
	fCaina_iz = Cafree_iz/(Cafree_iz + Kmact);	// Ca-dependency of the inactivation step

	// I1 inactivation gate (Na-dependent)
	alpha1_iz = (E1Na_iz*(fCaina_iz * a1on + (1.0 - fCaina_iz) * a1off));// *fixzero_slowVariable;	// rate from E1 to I1
	beta1_iz  = (fCaina_iz * b1on + (1.0 - fCaina_iz) * b1off); // *fixzero_slowVariable;	// rate from I1 to E1
	// I2 inactivation gate (Ca-dependent)
	alpha2_iz = (fCaina_iz*a2on + (1.0 - fCaina_iz)*a2off);	// rate from E1 to I2
	beta2_iz  = (fCaina_iz*b2on + (1.0 - fCaina_iz)*b2off);	// rate from I2 to E1
	alphaE_iz = k2*E2Na + k4*E2Ca;			// overall transition rate from E2 to E1
	betaE_iz  = k1*E1Na_iz + k3*E1Ca_iz;	// overall transition rate from E1 to E2

	E2NCX_iz = 1.0 - E1NCX_iz - I1NCX_iz - I2NCX_iz;
	INCX_iz = Frc_NCX_iz*maxINCX*( k1*E1Na_iz*E1NCX_iz - k2*E2Na*E2NCX_iz );
	INCXNa_iz = 3.0 * INCX_iz;
	INCXCa_iz = -2.0 * INCX_iz;


	// blk NCX
	E1Na_blk = 1.0 / (1.0 + pow((KmNai/Nai),3.0) * (1.0 + Cafree_blk/KmCai));	// Na binding probability
	E1Ca_blk = 1.0 / (1.0 + (KmCai/Cafree_blk) * (1.0 + pow((Nai/KmNai),3.0)));	// Ca binding probability
	fCaina_blk = Cafree_blk/(Cafree_blk + Kmact);	// Ca-dependency of the inactivation step

	// I1 inactivation gate (Na-dependent)
	alpha1_blk = (E1Na_blk*(fCaina_blk * a1on + (1.0 - fCaina_blk)*a1off));// *fixzero_slowVariable;// rate from E1 to I1
	beta1_blk  = (fCaina_blk*b1on + (1.0 - fCaina_blk)*b1off); // *fixzero_slowVariable;// rate from I1 to E1
	// I2 inactivation gate (Ca-dependent)
	alpha2_blk = (fCaina_blk*a2on + (1.0 - fCaina_blk)*a2off);	// rate from E1 to I2
	beta2_blk  = (fCaina_blk*b2on + (1.0 - fCaina_blk)*b2off);	// rate from I2 to E1
	alphaE_blk = k2*E2Na + k4*E2Ca;	// overall transition rate from E2 to E1
	betaE_blk  = k1*E1Na_blk + k3*E1Ca_blk;	// overall transition rate from E1 to E2

	E2NCX_blk = 1.0 - E1NCX_blk - I1NCX_blk - I2NCX_blk;
	INCX_blk = Frc_NCX_blk*maxINCX*( k1*E1Na_blk*E1NCX_blk - k2*E2Na*E2NCX_blk );
	INCXNa_blk = 3.0 * INCX_blk;
	INCXCa_blk = -2.0 * INCX_blk;

}

// Sodium-Potassium Pump

void comp_inak (double x[])
{
	//double frc = 1.0;
	//double scfmagNaK = 1.0;			// a scaling factor of Na/K pump
	//double dATP_NaK;
	//double tes1,tes2,tes3;
	double mult;
	double fVm,sfKdNaK;	// (VF/R/T)
	double KdNae,KdNai,KdKe,KdKi;	// Normalized Kd
	double Nai_c,Nae_c,Ki_c,Ke_c;		// normalized concentration referring to the dissociation constant
	double Nai_a,Nae_a,Ki_a,Ke_a;		// normalized concentration referring to the dissociation constant
	double a1p,a2p,a3p,a4p;
	double a1m,a2m,a3m,a4m; 			// transition rates in the reverse direction
	double denomi,numer;
	double VcycleC,VcycleA,Vcyc;
	double vm,Nai,Ki;
	vm    = x[0];
	Nai   = x[44];
	Ki    = x[45];

//	afNaK = baseNaK + deltaNaK * y_PKA;	// phospholemman of NaK phospholylated by PKA + base (activated by some other PKs) 
	//afNaK = 0.1 + 0.3*y_PKA;	// phospholemman of NaK phospholylated by PKA + base (activated by some other PKs) 

	sfKdNaK = 1.0;  // modulation by the phosphlemman
	//mult = -3/2;
	mult = -1.5;
	fVm   = vm*F/(R*T);					// VF/R/T
	KdNae = KdNae0*exp(dV_Nae*fVm);		// relative Kd
	KdNai = KdNai0*sfKdNaK*exp(dV_Nai*fVm);	// relative Kd for Nai
	KdKe  = KdKe0*pow(sfKdNaK,mult)*exp(dV_Ke*fVm);	// relative Kd for Ke
	KdKi  = KdKi0*exp(dV_Ki*fVm);		// relative Kd

	Nai_c = Nai / KdNai;		// relative concentration referring to the dissociation constant
	Nae_c = Nao / KdNae;		// relative concentration referring to the dissociation constant
	Ki_c  = Ki / KdKi;		// relative concentration referring to the dissociation constant
	Ke_c  = Ko / KdKe;		// relative concentration referring to the dissociation constant

	// transition rates in the forward direction
	a1p = (k1p * (Nai_c*Nai_c*Nai_c)) / (((1.0 + Nai_c)*(1.0+Nai_c)*(1.0+Nai_c)) + ((1.0 + Ki_c)*(1.0 + Ki_c)) - 1.0);
	a2p = k2p;
	a3p = k3p * (Ke_c*Ke_c) / (((1.0 + Nae_c)*(1.0 + Nae_c)*(1.0 + Nae_c)) + ((1.0 + Ke_c)*(1.0 + Ke_c)) - 1.0);
	a4p = k4p * MgATP_cyt / KdMgATP / (1.0 + MgATP_cyt / KdMgATP);
	
	// transition rates in the reverse direction
	a1m = k1m * MgADP_cyt;
	a2m = k2m * (Nae_c*Nae_c*Nae_c) / (((1.0 + Nae_c)*(1.0 + Nae_c)*(1.0 + Nae_c)) + ((1.0 + Ke_c)*(1.0 + Ke_c)) - 1.0);
	a3m = k3m * Pifree_cyt * H_cyt / (1.0 + MgATP_cyt / KdMgATP);
	a4m = k4m * (Ki_c*Ki_c) / (((1.0 + Nai_c)*(1.0 + Nai_c)*(1.0 + Nai_c)) + ((1.0 + Ki_c)*(1.0 + Ki_c)) - 1.0);

	denomi = (a1m + a1p) * a2m * a3m + a1p * a2p * (a3p + a3m) + a2p * a3p * (a4p + a4m) 
			+ (a2p + a2m) * a3m * a4m + (a1m + a1p) * a3p * a4p + a1m * (a3p + a3m) * a4m 
			+ a1p * (a2p + a2m) * a4p + a1m * a2m * (a4p + a4m);
	numer = a1p * a2p * a3p * a4p - a1m * a2m * a3m * a4m;

	VcycleC = numer / denomi;	// steady-state solution

	sfKdNaK = 0.72; // define the activated condition (Nai sensitivity is increased)
	KdNai = KdNai0*sfKdNaK*exp(dV_Nai*fVm);	// relative Kd for Nai
    KdKe = KdKe0*(pow(sfKdNaK,mult))*exp(dV_Ke*fVm);	// relative Kd for Ke
	//tes1=KdKe0;
	//tes2=pow(sfKdNaK,mult);
	//tes3=exp(dV_Ke*fVm);
	
	Nai_a = Nai / KdNai;		// relative concentration referring to the dissociation constant
	Nae_a = Nao / KdNae;		// relative concentration referring to the dissociation constant
	Ki_a  = Ki / KdKi;		// relative concentration referring to the dissociation constant
	Ke_a  = Ko / KdKe;		// relative concentration referring to the dissociation constant
//printf("Nai_a=%9.8lf,Nae_a=%9.8lf,Ki_a=%9.8lf,Ke_a=%9.8lf,KdKe=%9.8lf,aux=%9.8lf,aux2=%9.8lf,aux3=%9.8lf\n",Nai_a,Nae_a,Ki_a,Ke_a,KdKe,0.8*pow(sfKdNaK,mult)*exp(dV_Ke*fVm),exp(dV_Ke*fVm),pow(sfKdNaK,-1.5));
//printf("tes1=%9.8lf,tes2=%9.8lf,tes3=%9.8lf,tes4=%9.8lf\n",tes1,tes2,tes3,tes1*tes2*tes3);

	// transition rates in the forward direction
	a1p = (k1p * (Nai_a*Nai_a*Nai_a)) / (((1.0 + Nai_a)*(1.0+Nai_a)*(1.0+Nai_a)) + ((1.0 + Ki_a)*(1.0 + Ki_a)) - 1.0);
	a2p = k2p;
	a3p = k3p * (Ke_a*Ke_a) / (((1.0 + Nae_a)*(1.0 + Nae_a)*(1.0 + Nae_a)) + ((1.0 + Ke_a)*(1.0 + Ke_a)) - 1.0);
	a4p = k4p * MgATP_cyt / KdMgATP / (1.0 + MgATP_cyt / KdMgATP);
	
	// transition rates in the reverse direction
	a1m = k1m * MgADP_cyt;
	a2m = k2m * (Nae_a*Nae_a*Nae_a) / (((1.0 + Nae_a)*(1.0 + Nae_a)*(1.0 + Nae_a)) + ((1.0 + Ke_a)*(1.0 + Ke_a)) - 1.0);
	a3m = k3m * Pifree_cyt * H_cyt / (1.0 + MgATP_cyt / KdMgATP);
	a4m = k4m * (Ki_a*Ki_a) / (((1.0 + Nai_a)*(1.0 + Nai_a)*(1.0 + Nai_a)) + ((1.0 + Ki_a)*(1.0 + Ki_a)) - 1.0);

	denomi = (a1m + a1p) * a2m * a3m + a1p * a2p * (a3p + a3m) + a2p * a3p * (a4p + a4m) 
			+ (a2p + a2m) * a3m * a4m + (a1m + a1p) * a3p * a4p + a1m * (a3p + a3m) * a4m 
			+ a1p * (a2p + a2m) * a4p + a1m * a2m * (a4p + a4m);
	numer = a1p * a2p * a3p * a4p - a1m * a2m * a3m * a4m;

	VcycleA = numer / denomi;	// steady-state solution

	// summation of control and activated Vcycle
	Vcyc = afNaK * VcycleA + (1.0 - afNaK) * VcycleC;

  	INaK = maxINaK * Vcyc * rate_NaK;	// current  density (pA/pF)
	//INaK_Na = stoiNa * INaK;	// calculate ion flux
	//INaK_K = stoiK * INaK; // stoiNa = 3.0; stoiK = -2.0;
	//printf("INaK = %lf ",INaK);
	//printf("vcyc=%9.8lf,INaK=%9.8lf,INaK_Na=%9.8lf,INaK_K=%9.8lf\n",Vcyc,INaK,stoiNa*INaK,stoiK*INaK);

}

// Sarcolemmal Ca Pump,iz and blk (Ca-ATPase) 
void comp_ipmca (double x[])
{

	double Cafree_iz,Cafree_blk;
	Cafree_iz  = x[40];
	Cafree_blk = x[41];

	IPMCA_iz  = Frc_iz * maxIPMCA * pow(Cafree_iz,1.6) / (pow(KmPMCA,1.6) + pow(Cafree_iz,1.6));
	IPMCA_blk = Frc_blk* maxIPMCA * pow(Cafree_blk,1.6) / (pow(KmPMCA,1.6) + pow(Cafree_blk,1.6));

}

// K Background Current
void comp_ikb (double x[])
{
	//double scfIKb;
	//scfIKb = 1.0;

	IKbg = gkb*scfIKb*(x[0] - Ek);
}

// Na Background Current 

void comp_inab (double x[])
{
	INabg = gnab*(x[0] - Ena);
}

// Ca Background Current (iz and blk)
void comp_icab (double x[])
{
	ICabg_blk = Frc_blk*gcab*(x[0] - ECa_blk);
	ICabg_iz = Frc_iz*gcab*(x[0] - ECa_iz);
}

void comp_jup (double x[]) // SERCA
{
	double a1_pc,a2_pc,a3_pc;
	double a1_pa,a2_pa,a3_pa;
	double a1_mc,a2_mc,a3_mc;
	double a1_ma,a2_ma,a3_ma;
	double denomi_c,numer_c;
	double denomi_a,numer_a;
	double Vcycle_c,Vcycle_a;
	double JSERCA_control,JSERCA_active;
	double Kd_Cai,Kd_Cai_a,Kd_Casr;
	//double y_PKA,Cafree_blk,Ca_SRup;
	double Cafree_blk,Ca_SRup;

	//y_PKA      = x[49];
	Cafree_blk = x[41];
	Ca_SRup    = x[42];

//	afSERCA = baseSC + deltaSC * y_PKA;
//	afSERCA = 0.1 + 0.9*y_PKA;
// non-activated fraction of SERCA
	//orgKd_Cai = 0.91;	// original value Kd   huVEC 0.0027
	//orgKd_Casr = 2.24;	// original value Kd  huVEC   1.378

	//sfKdCa_SERCA = 1.0;	// control SERCA condition
	//Kd_Cai = orgKd_Cai * sfKdCa_SERCA
	Kd_Cai = 0.91 * 1.0;
	Kd_Casr = 2.24;

	//k1_plus = 25900	'mM-1 s-1    ***could be modified when the kinetic and dissociation constants are changed.
	//k2_plus = 2540	's-1
	//k3_plus = 20.5	's-1

	//k1_minus = 2	'mM-1 s-1     163.85               '
	//k2_minus = 67200	'mM-1 s-1
	//k3_minus = 149	'mM-1 s-1
	//a1_p = k1_plus*MgATP_cyt;
	//a2_p = k2_plus*pow(Cafree_blk/Kd_Cai,2)/(pow(Cafree_blk/Kd_Cai,2)*(1.0 + Hi*Hi/Kd_Hi) + Hi*Hi/ Kd_Hi*(1.0 + Hi/Kd_H1));
	//a3_p = k3_plus*Hsr*Hsr/Kd_Hsr/((Hi/Kd_H_release)*(1.0 + pow(Ca_SRup/Kd_Casr,2)) + (Hsr*Hsr/Kd_Hsr*(1.0 +Hi/Kd_H_release)));
	//a1_m = k1_minus*Hi*Hi/Kd_Hi/(pow(Cafree_blk/Kd_Cai,2)*(1.0 + Hi*Hi/Kd_Hi) + Hi*Hi/Kd_Hi*(1.0 + Hi/Kd_H1));
	//a2_m = k2_minus*MgADP_cyt*pow(Ca_SRup/Kd_Casr,2)*Hsr*Hsr/Kd_Hsr/((Hi/Kd_H_release)*(1.0 + (pow(Ca_SRup/Kd_Casr,2))) 
	//		+ (Hsr*Hsr/Kd_Hsr*(1.0 + Hi/Kd_H_release)));
	//a3_m = k3_minus * Pifree_cyt;
	a1_pc = 25900.0*MgATP_cyt;
	a2_pc = 2540.0*(Cafree_blk/Kd_Cai)*(Cafree_blk/Kd_Cai)/((Cafree_blk/Kd_Cai)*(Cafree_blk/Kd_Cai)*(1.0 + Hi*Hi/Kd_Hi) 
				+ Hi*Hi/ Kd_Hi*(1.0 + Hi/Kd_H1));
	a3_pc = 20.5*Hsr*Hsr/Kd_Hsr/((Hi/Kd_H_release)*(1.0 + (Ca_SRup/Kd_Casr)*(Ca_SRup/Kd_Casr)) 
				+ (Hsr*Hsr/Kd_Hsr*(1.0 +Hi/Kd_H_release)));
	a1_mc = 2.0*Hi*Hi/Kd_Hi/((Cafree_blk/Kd_Cai)*(Cafree_blk/Kd_Cai)*(1.0 + Hi*Hi/Kd_Hi) + Hi*Hi/Kd_Hi*(1.0 + Hi/Kd_H1));
	a2_mc = 67200.0*MgADP_cyt*(Ca_SRup/Kd_Casr)*(Ca_SRup/Kd_Casr)*Hsr*Hsr/Kd_Hsr/
			((Hi/Kd_H_release)*(1.0 + (Ca_SRup/Kd_Casr)*(Ca_SRup/Kd_Casr)) + (Hsr*Hsr/Kd_Hsr*(1.0 + Hi/Kd_H_release)));
	a3_mc = 149.0*Pifree_cyt;

	denomi_c = a2_pc*a3_pc + a1_mc*a3_pc + a1_mc*a2_mc + a1_pc*a3_pc + a2_mc*a1_pc 
								+ a2_mc*a3_mc + a1_pc*a2_pc + a3_mc*a1_mc + a3_mc*a2_pc;
	numer_c  = a1_pc*a2_pc*a3_pc - a1_mc*a2_mc*a3_mc;
	Vcycle_c = numer_c/ denomi_c;
	
	JSERCA_control = ((1.0 - afSERCA)*maxISERCA*Vcycle_c)/(2.0*F)*freezeSR;	// betafactor
	//dATP_SERCA_control = JSERCA_control/2.0; 			// stoichiometry = 2Ca/1ATP      amole / ms
	
// activated fraction of SERCA
	
	//sfKdCa_SERCA = 0.5;	// activated condition by the beta stimulation
	//Kd_Cai = orgKd_Cai * sfKdCa_SERCA
	Kd_Cai_a = 0.91*0.5;
	// Kd_Casr = 2.24;
	//k1_minus = 2/pow(sfKdCa_SERCA, 2);
	//k1_minus = 2/0.5/0.5;
	a1_pa = 25900.0*MgATP_cyt;
	a2_pa = 2540.0*pow(Cafree_blk/Kd_Cai_a,2)/(pow(Cafree_blk/Kd_Cai_a,2)*(1.0 + Hi*Hi/Kd_Hi) + Hi*Hi/ Kd_Hi*(1.0 + Hi/Kd_H1));
	a3_pa = 20.5*Hsr*Hsr/Kd_Hsr/((Hi/Kd_H_release)*(1.0 + (pow(Ca_SRup/Kd_Casr,2))) + (Hsr*Hsr/Kd_Hsr*(1.0 +Hi/Kd_H_release)));
	a1_ma = (2.0/0.5/0.5)*Hi*Hi/Kd_Hi/(pow(Cafree_blk/Kd_Cai_a,2)*(1.0 + Hi*Hi/Kd_Hi) + Hi*Hi/Kd_Hi*(1.0 + Hi/Kd_H1));
	a2_ma = 67200.0*MgADP_cyt*pow(Ca_SRup/Kd_Casr,2)*Hsr*Hsr/Kd_Hsr/((Hi/Kd_H_release)*(1.0 + (pow(Ca_SRup/Kd_Casr,2))) 
			+ (Hsr*Hsr/Kd_Hsr*(1.0 + Hi/Kd_H_release)));
	a3_ma = 149.0*Pifree_cyt;
	
	denomi_a = a2_pa*a3_pa + a1_ma*a3_pa + a1_ma*a2_ma + a1_pa*a3_pa + a2_ma*a1_pa
								+ a2_ma*a3_ma + a1_pa*a2_pa + a3_ma*a1_ma + a3_ma*a2_pa;
	numer_a  = a1_pa*a2_pa*a3_pa - a1_ma*a2_ma*a3_ma;
	Vcycle_a = numer_a / denomi_a;
	
	JSERCA_active = (afSERCA*maxISERCA*Vcycle_a)/(2.0*F)*freezeSR; 	// betafactor 
	//dATP_SERCA_active = JSERCA_active/2.0;			// stoichiometry = 2Ca / 1ATP      amole / ms

// summation of control and activated fractions
	JSERCA = JSERCA_control + JSERCA_active;	// sum of active and non-activated components
	//dATP_SERCA = dATP_SERCA_control + dATP_SERCA_active;	// ATP consumption rate determined with the stoichiometry = 2Ca / 1ATP      amole / ms

}

void comp_jtr (double x[])
{
	double Cafree_SRrl,Ca_SRup;

	Ca_SRup     = x[42];
	Cafree_SRrl = x[43];

	Jtrans_SR = Ptrans * (Ca_SRup - Cafree_SRrl);

}

void comp_CaRU (double x[])
{
	double vm,fVm_Act,fVm_inAct,mEtan12,expF;
	double pC1,pC2;
	double po_LCC,po_RyR,pot_RyR;
	double Ca_ndLR,Ca_nd0R,Ca_ndL0,Ca_nd00;
	double kco,koc,Akco,kRco1,kRoc1;
	//double kco2,koc2,Akco2,kRco2,kRoc2;
	double kco2,   Akco2,kRco2,kRoc2;
	double ft,ft_rest,kco_rest;
	double ft2,ft_rest2,kco_rest2;
	//double kvco,kvoc;
	//double kcaco,kcaoc;
	double perNa_ICaL,perK_ICaL;
	//double Yooo,Yooc,Ycoo,Ycoc,Ycco,Yoco,Yocc;
	double Yooo,Yooc,Ycoo,Ycco,Yoco;
	double Cafree_jnc,Cafree_iz,Cafree_SRrl;
	double tPCaL;

	         vm = x[0];
	       Yooo = x[27];
	       Yooc = x[28];
	       Ycoo = x[29];
	       Ycco = x[31];
	       Yoco = x[32];
	 Cafree_jnc = x[39];
	 Cafree_iz  = x[40];
	Cafree_SRrl = x[43];

	//fVmAct    = 1.0/(3.734*exp(-(vm + VshiftLCCact)/8.5) + 0.35*exp(-(vm + VshiftLCCact)/3500.0));
	//fVmInAct  = 1.0/(4.65*exp((vm + VshiftLCCinact)/15.0) + 1.363*exp((vm + VshiftLCCinact)/100.0));
	fVm_Act    = 1.0/(3.734*exp(-vm/8.5) + 0.35*exp(-vm/3500.0));
	fVm_inAct  = 1.0/(4.65*exp(vm/15.0) + 1.363*exp(vm/100.0));
	mEtan12   = 1.0/(8084.0*exp(vm/10.0) + 158.0*exp(vm/1000.0)) + 1.0/(134736.0*exp(-vm/5.0) + 337.0*exp(-vm/2000.0));

	po_LCC = Yooo + Yooc;				// open probability of LCC
	po_RyR = Yooo + Ycoo + Ycco + Yoco;	// open probability of RyR
	
	expF   = exp(-vm/RTon2F);
	Ca_ndLR = (Cafree_jnc + Cafree_SRrl*JRyR/gD_nd + expF*Cao*JLCC/gD_nd*(vm/RTon2F)/(1.0 - expF))/(1.0 + JRyR/gD_nd + JLCC/gD_nd*vm/RTon2F/(1.0 - expF));
	Ca_nd0R = (Cafree_jnc + Cafree_SRrl*JRyR/gD_nd)/(1.0 + JRyR/gD_nd);
	Ca_ndL0 = (Cafree_jnc + expF*Cao*JLCC/gD_nd*(vm/RTon2F)/(1.0 - expF))/(1.0 + JLCC/gD_nd*vm/RTon2F/(1.0 - expF));	
	Ca_nd00 = Cafree_jnc;
	Ca_Noise = Ca_nd00; // Ca_nd00*rand(): from protocolGeneral2.vb

/*	if(fabs(vm)<1E-4){
		Ca_ndLR = (Cafree_jnc + Cafree_SRrl*JRyR/gD_nd + Cao*JLCC/gD_nd)/(1.0 + JRyR/gD_nd + JLCC/gD_nd);
	} else {
		Ca_ndLR = (Cafree_jnc + Cafree_SRrl*JRyR/gD_nd + expF*Cao*JLCC/gD_nd*(vm/RTon2F)/(1.0 - expF))/(1.0 + JRyR/gD_nd + JLCC/gD_nd*vm/RTon2F/(1.0 - expF));
	}
*/
	// Ca_nd with L open
/*	if(fabs(vm)<1E-4){
		Ca_ndL0 = (Cafree_jnc + Cao*JLCC/gD_nd)/(1.0 + JLCC/gD_nd);
	} else {
		Ca_ndL0 = (Cafree_jnc + expF*Cao*JLCC/gD_nd*(vm/RTon2F)/(1.0 - expF))/(1.0 + JLCC/gD_nd*vm/RTon2F/(1.0 - expF));	
	}
*/

	kco = 3.0*0.4/(1.0 + pow(0.025/Ca_ndL0,2.7));
	koc = 3.0*0.5564;
	ft = kco/(kco + koc);

	Akco = 3.0*0.4*(0.1 + Cafree_SRrl)/(1.0 + pow(0.025/Ca_ndLR,2.7));
	kRco1 = 7.0*(ft * Akco);

	kco_rest = 3.0*0.4/(1.0 + pow(0.025/Ca_nd00,2.7));
	ft_rest  = kco_rest/(kco_rest + koc);
	pC1 = koc/(ft_rest*Akco + koc);
	kRoc1 = koc*pow(pC1,(9.0*0.74));

	kYoooYooc = kRoc1;	// trigger step closing rate
	kYoocYooo = kRco1;	// trigger step opening rate

	// ******  Case of LCC is closed ********
	kco2 = 3.0* 0.4/ (1.0 + pow(0.025/Ca_Noise,2.7));
	ft2 = kco2/(kco2 + koc);
	Akco2 = 3.0*0.4*(0.1 + Cafree_SRrl)/(1.0 + pow(0.025/Ca_nd0R,2.7));
	kRco2 = 7.0*(ft2 * Akco2);
	kco_rest2 = 3.0*0.4/(1.0 + pow(0.025/Ca_nd00,2.7));
	ft_rest2  = kco_rest2/(kco_rest2 + koc);
	pC2 = koc/(ft_rest2*Akco2 + koc);
	kRoc2 = koc*pow(pC2,(9.0*0.74));

	kYcooYcoc = kRoc2;	// trigger step closing rate     Ycoo
	kYcocYcoo = kRco2;	// trigger step opening rate
	
	kYccoYccc = kRoc2;	// trigger step closing rate   Ycco
	kYcccYcco = kRco2;	// trigger step opening rate
	
	kYocoYocc = kRoc2;	// trigger step closing rate   Yoco
	kYoccYoco = kRco2;	

	kYcooYooo = fVm_Act;
	kYoooYcoo = fVm_inAct;
	kYcocYooc = fVm_Act;
	kYoocYcoc = fVm_inAct;

	kYccoYoco = fVm_Act;
	kYocoYcco = fVm_inAct;
	kYcccYocc = fVm_Act;
	kYoccYccc = fVm_inAct;

	kYoooYoco = (Ca_ndLR/KL)*fVm_Act/TL;
	kYocoYooo = mEtan12;

	kYoocYocc = (Ca_ndL0/KL)*fVm_Act/TL;
	kYoccYooc = mEtan12;
	
	kYcooYcco = (Ca_nd0R/KL)*fVm_Act/TL;
	kYccoYcoo = mEtan12;

	kYcocYccc = (Ca_nd00/KL)*fVm_Act/TL;
	kYcccYcoc = mEtan12;
	
	//afLCC2 = 1.0 + 1.3*y_PKA;
	tPCaL = (oPCaL * scf_LCC) * afLCC;

	GHK_Ca_LR  = zca*vm/RTonF*(Ca_ndLR - Cao*exp(-zca*vm/RTonF))/(1.0 - exp(-zca*vm/RTonF));
	GHK_Ca_L0  = zca*vm/RTonF*(Ca_ndL0 - Cao*exp(-zca*vm/RTonF))/(1.0 - exp(-zca*vm/RTonF));
	ICaLCa_LR = Frc_CaL_jnc * tPCaL * GHK_Ca_LR * Yooo * ATPfactor;
	ICaLCa_L0 = Frc_CaL_jnc * tPCaL * GHK_Ca_L0 * Yooc * ATPfactor;
	//printf("ICaLCa_LR=%lf,Frc=%lf,tPCaL=%lf,GHK_Ca_LR=%lf,Yooo=%lf,ATP=%lf\n",ICaLCa_LR,Frc_CaL_jnc,tPCaL,GHK_Ca_LR,Yooo,ATPfactor);
	perNa_ICaL = 0.0000185 * tPCaL;
	perK_ICaL = 0.000367 * tPCaL;
	
	ICaLNa_jnc = Frc_CaL_jnc * perNa_ICaL * GHK_Na * po_LCC * ATPfactor;
	ICaLK_jnc  = Frc_CaL_jnc * perK_ICaL  * GHK_K  * po_LCC * ATPfactor;
	
	pot_RyR = (1.0 - 0.00006)*po_RyR + 0.00006;	
	
	Jrel_SR = PRyR*pot_RyR*(Cafree_SRrl - Cafree_jnc);

	Jleak_SR = PRyR * 0.00006*(Cafree_SRrl - Cafree_iz);

}

void comp_Ca_buff_jnc (double x[]) // x[42]
{
	double f1,f2;
	double Cafree_jnc;
	Cafree_jnc = x[39];

	f1 = KdL_iz*BtotL_jnc / (Cafree_jnc + KdL_iz) / (Cafree_jnc + KdL_iz);
	f2 = KdH_iz*BtotH_jnc / (Cafree_jnc + KdH_iz) / (Cafree_jnc + KdH_iz);

	Buf_jnc = 1.0/(1.0 + f1 + f2);
}

void comp_Ca_buff_iz (double x[]) // x[43]
{
	double f1,f2;
	double Cafree_iz;
	Cafree_iz = x[40];

	f1 = KdL_iz*BtotL_iz / (Cafree_iz + KdL_iz) / (Cafree_iz + KdL_iz);
	f2 = KdH_iz*BtotH_iz / (Cafree_iz + KdH_iz) / (Cafree_iz + KdH_iz);

	Buf_iz = 1.0/(1.0 + f1 + f2);

}

void comp_Ca_buff_blk (double x[]) // x[44]
{
	double Cafree_blk;
	double f1,f2,f3;

	Cafree_blk = x[41];

	f1 = BtotCaM*KoffCaM/KonCaM/ (Cafree_blk + KoffCaM/KonCaM) / (Cafree_blk + KoffCaM/KonCaM);
	f2 = BtotTnCh*KoffTnChCa/KonTnChCa / (Cafree_blk + KoffTnChCa/KonTnChCa) / (Cafree_blk + KoffTnChCa/KonTnChCa);
	f3 = BtotSR*KoffSR/KonSR / (Cafree_blk + KoffSR/KonSR) / (Cafree_blk + KoffSR/KonSR);

	Buf_blk = 1.0/(1.0 + f1 + f2 + f3);
}

void comp_buffSR_rl (double x[]) // x[**]
{
	double Cafree_SRrl;
	double ff;

	Cafree_SRrl = x[43];

	ff = KdCsqnCa*Btot_Csqn / (Cafree_SRrl + KdCsqnCa) / (Cafree_SRrl + KdCsqnCa);

	Buf_SRrl = 1.0/(1.0 + ff);
}

void comp_bund_dif (double x[])
{
	double Cafree_jnc,Cafree_iz,Cafree_blk;

	Cafree_jnc = x[39];
	Cafree_iz  = x[40];
	Cafree_blk = x[41];

	JCa_jnciz = GCa_jnciz*(Cafree_jnc - Cafree_iz);
	JCa_izblk = GCa_izblk*(Cafree_iz  - Cafree_blk);

}

void comp_IP3R (double x[])
{
	double Cafree_jnc,Cafree_SRrl;
	double O_ip3,A_ip3;
	double pO_IP3R,uMCa;

	Cafree_jnc  = x[39];
	Cafree_SRrl = x[43];
	 O_ip3 = x[35];
	 A_ip3 = x[38];

	uMCa = 1E+3*Cafree_jnc;	// not original source
//	uMCa = 1E-3*Cafree_jnc;	//Ca concentration microM

	phi_1 = (k_1a * L_1 + l_2a) * uMCa / (L_1 + uMCa * (1.0 + L_1 / L_3));
	phi_2 = (k_2a * L_3 + l_4a * uMCa) / (L_3 + uMCa * (1.0 + L_3 / L_1));
	phi_2b = (k_2b + l_4b * uMCa) / (1.0 + uMCa / L_5);
	phi_3 = k_3a * L_5 / (uMCa + L_5);
	phi_4 = (k_4a * L_5 + l_6a) * uMCa / (uMCa + L_5);
	phi_4b = L_1 * (k_4b + l_6b) / (uMCa + L_1);
	phi_5 = (k_1a * L_1 + l_2a) * uMCa / (uMCa + L_1);
	
	pO_IP3R = pow((0.1*O_ip3 + 0.9*A_ip3),4.0);		// calculate pO_IP3R in response to [IP3]

	J_ip3R = Pip3_wholeCell*pO_IP3R*(Cafree_SRrl - Cafree_jnc)*freezeSR;
}

void comp_aAR_kinetics (double x[])
{
//'****************** alpha AR kinetics ***************
	double conIP3_cyt;
	
	//conIP3_cyt = 0.015;

	//conIP3 = conIP3_cyt;	//'concentration of IP3 at a given time of integration
	//if (conIP3 < 0.015)  conIP3 = 0.015;	// 'lower limit of [IP3] in the cytosol

//'***************** IKb activation ************************
	//scfmagIKb = scfmagIKb + (scfmagIKb_inf - scfmagIKb) * (dt / tauIKbg);     'time-dependent activation of  IKbg

}

void comp_bAR_kinetics (double x[])
{
	double ARS464,ARS301;
	double GsaGTP;
  //double Gsbg;
	double GsaGDP;
	double cAMP;
	double cataliticC;
	double y_PKA;
	double GsaGTPtot;

	ARS464    = x[46]; // desensitize rate by phosphorylation by bBarK
	ARS301    = x[47]; // desensitize rate by PKA
	GsaGTP    = x[48]; // time course of total G-GTP complex
	//Gsbg    = x[49]; // time course of Gs betagamma- complex
	GsaGDP    = x[50]; // time course of Gs alpha-GDP complex,  (hydrolysis of GTP)
	cAMP      = x[51]; // free cyclic AMP concentration in cytosol
	cataliticC= x[52]; // Cat (Cat4) in Cat45
	y_PKA     = x[53];

	// active beta adrenoceptor after subtraction of  desensitized conformations
	// ARtot: total b1-adrenergic receptor, 0.0000132 (mM)
	ARact = ARtot - (ARS464 + ARS301);

	// Gs, GTP binding protein, Gsfree + GTP-bound Gs + GDP-bound Gs
	// Gstot: total Gs protein, 0.00383 (mM)
	//Gsfree = Gstot - GsaGTPtot - GsaGDP;
	GsaGTPtot = GsaGTP + ACtot*GsaGTP/(GsaGTP + KdACG);
	Gsfree = Gstot - GsaGTPtot - GsaGDP;
	
	B2 = 1.0/(1.0 + ACtot*KdACG/(GsaGTP + KdACG)/(GsaGTP + KdACG));

	// concentration of GsGTP-AC complex
	//AC_GaGTP = GsaGTPtot - GsaGTP;
	AC_GaGTP = ACtot*GsaGTP/(GsaGTP + KdACG);
	// concentration of free AC
	AC = ACtot - AC_GaGTP;
	
	// cAMP-activated protein kinase   Cat123 and Cat45
	// PKAtot = 0.001; total protein kinase A (dimer), 0.001 (mM)
	// Cattot = PKAtot*2.0 = Cat123 + Cat45;
	// Cat123 = Cattot - Cat45;
	Cat45 = cataliticC + PKItot*cataliticC/(Kd4 + cataliticC);
	Cat123 = 2.0*0.001 - Cat45;
	C_PKI = PKItot*cataliticC/(Kd4 + cataliticC); // PKA - PKI complex (mM)

	// [cAMP] translate the cAMP_free function to a state [cAMP]free
	//A2R = Cat45;
	//cAMP = cAMPtot - ARC_bAR - 2.0*A2RC - 2.0*A2R;
	B10 = -Kd12*Cat123*cAMP;
	B11 = +Kd12*Cat123*cAMP*cAMP*(2.0+2.0*Kd3/cataliticC);
	B12 = 2.0*Kd12*Kd12*Cat123*cAMP*(2.0+2.0*Kd3/cataliticC);
	B13 = Kd12*Kd12*Kd12*Cat123;
	B1  = 1.0/(1.0 + (B10 + B11 + B12 + B13)/(cAMP*cAMP + Kd12*cAMP + Kd12*Kd12));

	// instantaneous distribution of AR states after the application of a Ligand 
	ARfree = ARact / (1.0 + Ligand / KdLR + Gsfree * Ligand / KdGLR / KdLR + Gsfree / KdGR);
	LAR = Ligand * ARfree / KdLR;
	LARGs = ARfree * Gsfree * Ligand / KdGLR / KdLR;
	ARGs = ARfree * Gsfree / KdGR;

	// PKA (activation)
	// alpha and beta for the original Kd = 0.009 of the instantaneous equilibrium A2RC + Cat = CA2RC
 	beta = alpha * Kd3;

	// A2RC in Cat123
	A2RC = Cat123/((Kd12/cAMP)*(Kd12/cAMP)+(Kd12/cAMP)+1.0);
	ARC_bAR = A2RC*Kd12/cAMP;
	RC = ARC_bAR*Kd12/cAMP;

	alpha_PKA = alpha * cataliticC; // = beta * A2R * Cat / Cat45     where A2R = Cat45 (x[48]) (ms-1)
	beta_PKA = beta/( (Kd12/cAMP)*(Kd12/cAMP) + (Kd12/cAMP) + 1.0); // ms-1  = alpha * A2RC

	B0 = 1.0/(1.0 + Kd4*PKItot/(Kd4 + cataliticC)/(Kd4 + cataliticC));

	// Action of the catalytic subunit, Cat
	// under an assumption that the proportion of phosphorylated target protein (y_PKA) changes with a tau of 40 sec.
	// for example, the spontaneous rate of APDAD evoked by NA application takes an exponential time course 
	// with a time constant of approximately 40 sec.

	kFcat = Kf * cataliticC;
	// kBpp = Kb * [PP]      note 1/tau = alpha + beta
	kBpp = 1.0/tau_phosphorylation - kFcat;

	// time-dependent active fraction of target protein phosppholylated *********************
	// LCC phosppholylated by PKA + base (activated by some other PKs)
	afLCC = baseCaL + deltaCaL * y_PKA;
	// phospholamban of SERCA phospholylated by PKA + base (activated by some other PKs)
	afSERCA = baseSC + deltaSC * y_PKA;
	// phospholemman of NaK phospholylated by PKA + base (activated by some other PKs)
	afNaK = baseNaK + deltaNaK * y_PKA;
	
}
