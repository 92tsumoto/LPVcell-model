
#include "syspara.h"

void val_consts(FILE *fp1)
{
	//int i,w;
	//double v_old,dvdt,dvdt_new;

	// frequency
	var.omega = 1.0/var.BCL;
	//var.omega = 2.0*M_PI/var.BCL;

	// Cell Geometry */
	// Forbes M, Sperelakis N. Ultrastructure of mammalian cardiac muscle. 
	// In: Sperelakis N, ed. Physiology and pathophysiology of the heart.
	// 2nd edition. Boston, MA: Kluwer Academic; 1989:3-41.
	Sc     = 0.59;			// scaling factor inducing from HuVEC model
	Vcell  = 22373.0;       // Cell Volume (pL): 120x37.62x8.4xscaling-factor
	Vblk   = Vcell*0.68;    // Bulk(myoplasm) volume (pL)
	Viz    = Vcell*0.035;    // IZ volume (pL)
	Vjnc   = Vcell*0.008;    // JSR volume (pL)
	Vcyt   = Vjnc + Viz + Vblk;	// cytosolic space (pF)
	Vsr    = Vcell*0.072;    // whole SR volume (pL) 
	Vsr_rl = Vsr*0.225;      // fractional volume of the junctional SR (pL)
	Vsr_up = Vsr*0.775;      // fractional volume of the network SR (pL)
	//var.ageo = 2*M_PI*var.a*var.a+2.0*M_PI*var.a*var.length;  // eometric membrane area: 7.671e-5 (cm^2)
	//var.acap = var.ageo*2;			// Capacitive membrane area: 1.534e-4 cm^2 (cm^2)
	Cm = 192.46*Sc;						// 113.5514 (pF):Capacitive membrane area (pF)

// only for the mitochondria model
//	Rcm_mit= 4.35;				// volume ratio of cytosol/mitochondria (Vol_cyt / Vol_mit):  
//	Vmit   = Vcyt/Rcm_mit;		// Mitochondria volume (pL)
//	Cm_mit = 0.001*Vmit*F;		// Mitochondria surface membrane capacitance (pF)
	//ampRatio = Vcyt*Sc;			// Vcyt * vRatio_mit <-- ? From Umehara et al., IJMS, 2019
//	ampRatio = 110036.0*Sc;			// Vcyt * vRatio_mit <-- ? From Umehara et al., IJMS, 2019

// Fraction of Currents jnc, iz, blk
	Frc_iz = 0.1;				// fractional distribution at the iz myoplasm for ion transporters except LCC and NCX
	Frc_blk = 0.9;				// fractional distribution at the blk myoplasm for ion transporters except LCC and NCX
	Frc_CaL_jnc = 0.75;			// fractional distribution of ICaL at the junctional myoplasm
	Frc_CaL_iz  = 0.15;			// fractional distribution of ICaL at the iz myoplasm
	Frc_CaL_blk = 0.10;			// fractional distribution of ICaL at the blk myoplasm
	Frc_NCX_jnc = 0.03;			// fractional distribution of NCX at the junctional myoplasm
	Frc_NCX_iz  = 0.25;			// fractional distribution of NCX at the iz myoplasm
	Frc_NCX_blk = 0.72;			// fractional distribution of NCX at the blk myoplasm

// Ion Valences
	zna = 1;  // Na valence
	zk = 1;   // K valence
	zca = 2;  // Ca valence

// invariant constant
	RTonF  = R*T/F;
	RTon2F = R*T/(zca*F);

// Extracellular Concentrations
	Nao   = 140.0;	// Initial Bulk Medium Na (mM)
	Ko    = 5.4;	// Initial Bulk Medium K (mM)
	Cao   = 1.8;	// Initial Bulk Medium Ca (mM)

// Cytosol Pi, AMP, ADP free, MgADP, ATP free, MgATP, Cr,
	// pH_cyt = -Log10(H_cyt / 1000)
	// HBuff_cyt = 25 / 2.3 / H_cyt
	// Pifree_cyt = Pi_cyt / (1 + Math.Pow(10, (pH_cyt - pKa_Pi)))
	// AMP_cyt = AXPSUM_cyt - (ATPt_cyt + ADPt_cyt), mM, AMP concentration  total
	// ADPfr_cyt = ADPt_cyt / (1 + Mg_cyt / KdADP_cyt), mM, Instantaneous equilibrium for Mg + ADP binding
	// MgADP_cyt = ADPt_cyt - ADPfr_cyt
	// ATPfr_cyt = ATPt_cyt / (1 + Mg_cyt / KdATP_cyt), mM,  Instantaneous equilibrium for Mg + ATP binding
	// MgATP_cyt = ATPt_cyt - ATPfr_cyt
	// Cr_cyt = CrSUM_cyt - PCr_cyt, mM, creatine concentration
	Pifree_cyt = 0.50872066859173026;  
	AMP_cyt = 0.00033459021041526427;
	ADPfr_cyt = 0.0022536111580301241;
	MgADP_cyt = 0.025978226605534577;
	ATPfr_cyt = 0.039789862258604494;
	MgATP_cyt = 6.631643709767415;
	Cr_cyt = 12.6772372798697;
	H_cyt  = 0.0001; // (mM) cytosol H concentration   Constant at present

// Fast sodium current
	kI2O = 0.0001312;
	//kI1I2 = 0.00534;
	//kI1O = 0.01;
	//fL = 0.13125;
	//PNa = 8.1072; // Himeno et al., BJ, 2015
	PNa = 9.4584; // Himeno et al., BJ, 2015
	rPK_Ina = 0.15;	// Umehara et al., IJMS, 2019
	rate_INa = 1.0; // Umehara et al., IJMS, 2019

// IKur current // (nS/pF).
	GKur = 0.0672;     // 0.16*0.6*0.7 ?   G_Kur = 0.0672 nS by Umehara et al., IJMS,2019	
	//rate_Kur = 1.0/Cm;
	rate_Kur = 1.0;

// Transient outward current (mS/uF)
	GKto = 0.19318*0.7;  // '0.02975*0.45*1000/Cm  'microS {init: 0.02975};
	fast_ito = 0.0;		// fraction of fast component Ito
	slow_ito = 1.0;		// fraction of slow component Ito
	rate_Ito = 1.0;

// L-type calcium current
	KL = 0.00154;	//0.00154 = 4.4*0.35/1000, (mM): Ca-sensivity is increased by scaling 0.35
	TL = 147.51;
//	gD_free = 0.065;	// (um3ms-1) Diffusion const between channel mouth and  iz-space or blk_space
//	JLCC = 0.000913;	// (um3ms-1) conductance of a unit LCC

	ATPCaL = 6.0;		// a concentration of ATP for the CaL gating
	ATPfactor = 1.0/(1.0+((1.4/ATPCaL)*(1.4/ATPCaL)*(1.4/ATPCaL)));
	oPCaL  = 14.21*0.42*1.5*1;	// about 60% of Umehara (original 14.21)
	//tPCaL; // time-dependent PCaL, about 60% of Umehara (original 14.21)
	scf_LCC = 1.0; // scaling factor for the number of whole cell ICaL channel
	
	VshiftLCCact = 0.0;
	VshiftLCCinact = 0.0;

// time-dependent active fraction of target protein phosppholylated 
	baseCaL  = 1.0;	// basal fraction of active (phosphorylated) LCC
	deltaCaL = 1.3;	// non-phosphorylated faction of LCC

	baseSC   = 0.1; // basal fraction of active (phosphorylated) phospholamban (SERCA)
    deltaSC  = 0.9; // non-phosphorylated faction of  phospholamban (SERCA)

	baseNaK  = 0.1;	// basal fraction of activated (phosphorylated) phospholemman (Na/K)
	deltaNaK = 0.3;	// non-phosphorylated faction of phospholemman (Na/K)

// Rapid delayed rectifier potassium current (Ikr) (mS/uF)
	gKr = 0.06839*4.0;
	qKr = pow((Ko/4.5),0.2);
	rate_Kr = 1.0;

// Inward rectifier K current: Ik1 (mS/uF)
	frac_mode_K1 = 0.9;
	SPM = 0.005*1000.0;		// SPM should be given in micro M in this model
	GK1 = 0.04546*0.7;                //revised by AN on 190405
	qK1 = pow((Ko/4.5),0.4);
	fK1 = qK1/(1.0 + exp(-(Ko - 2.2)/0.6));
	
// Cloride ion currnet I_Clh: Okamoto et al., JMCC, 2012,2014
  //GClh = 0.06*0.0;
	GClh = 0.06;
	rate_Clh = 1.0;

// Sodium-Calcium Exchanger (NCX)
	rate_NCX_jnc = 1.0;
	rate_NCX_iz = 1.0;
	rate_NCX_blk = 1.0;
	maxINCX = 10.99; 	// 61.06*0.5*0.3*1.2, reduced by 0.18
	a1off = 0.0015;
	a1on = 0.002;
	b1off = 0.0000005;
	b1on = 0.0012;
	a2off = 0.02;
	a2on = 0.00006;
	b2off = 0.0002;
	b2on = 0.18;
	// parameters for the Ca- and Na-binding to transport sites
	KmNao = 87.5; // (mM)
	KmNai = 20.74854;
	KmCao = 1.38;
	KmCai = 0.0184;
	Kmact = 0.004;
	//nHNa = 3.0;
	//partition = 0.32;
	k3 = 1.0;
	k4 = 1.0;
	E2Na = 1.0 / (1.0 + pow((KmNao / Nao),3.0) * (1.0 + Cao / KmCao));  // Na binding probability
	E2Ca = 1.0 / (1.0 + (KmCao / Cao) * (1.0 + pow((Nao / KmNao),3.0)));    // Ca binding probability

// Na-K Pump
	k1p = 0.72;
	k1m = 0.08;
	k2p = 0.08;
	k2m = 0.008;
	k3p = 4.0;
	k3m = 8000.0;
	k4p = 0.3;
	k4m = 0.2;

	KdKe0   = 0.8;
	KdKi0   = 18.8;
	KdNae0  = 26.8;
	KdNai0  = 5.0;
	KdMgATP = 0.6;

	stoiNa = 3.0;
	stoiK  = -2.0;

	dV_Nae = 0.44;
	dV_Nai = -0.14;
	dV_Ke  = 0.23;
	dV_Ki  = -0.14;
	//sfKdNaK = 1;	// modulation by the phosphlemman
	//************activated fraction of Na/K pump**************
	//sfKdNaK = 0.72;	// define the activated condition (Nai sensitivity is increased)

	maxINaK = 36.382;	// 25.1779*1.7*0.85, increased by x1.445
	rate_NaK = 1.0;		// a scaling factor of Na/K pump

// Sarcolemmal Ca Pump (PMCA)
	maxIPMCA = 0.19;		// Max. Ca current through sarcolemmal Ca pump (mS/uF)
	KmPMCA = 0.0005;		// Half-saturation concentration of sarcolemmal Ca pump (mM)

// K Background Current (mS/uF) 
	gkb = 0.025883;
	tauIKbg = 80000.0;		//'30 sec activation of IKbg during the alpha stimulation
	//scfIKb_inf = 1.0;

// Ca Background Current 
	gcab = 0.00014028;			// nS 

// Na Background Current 
	gnab = 0.0037206;

// SR calcium release flux, via RyR (Jrel)
	JRyR    = 0.02;		// (um3ms-1) conductance of a unit RyR in a single CaRU 
	JLCC    = 0.000913;	// (um3ms-1) conductance of a unit LCC
	gD_nd   = 0.065;	// (um3ms-1) Diffusion const between channel mouth and nano domain_space
	gD_free = 0.065;	// (um3ms-1) Diffusion const between channel mouth and iz-space or blk_space

// calcium uptake via SERCA pump (Jup)
	maxISERCA = 1341.204*Sc;// 447068*3*Sc_Cell, 106444.8*Sc_Cell*0.7*6 
	Kd_H1 = 1.09E-5;	// (mM) binding H to Ei , inhibition of Ca binding
	Kd_Hi  = 3.54E-3;	// (mM^2) bomdoing [H]i to transport
	Kd_Hsr = 1.05E-8;	// (mM^2) binding [H]sr  to transport
	Kd_H_release = 7.24E-5;	// (mM) release of H from the SERCA
	Hi  = 1E-4;
	Hsr = 1E-4;
	freezeSR = 1.0;

// Translocation of Ca Ions from NSR to JSR
	Ptrans = 4.8037*Sc;		// 22.745*0.1*Sc_Vol:
							// The transfer rate is much reduced to ensure the slow recovery of
							// Ca content in the junctional SR space.
	PRyR   = 4177.369*Sc;	// 5967.67 * Sc_Cell * 0.7

// aAR 
	m_inf = 1.0;
	h_inf = 0.0;
	tau_mi = 10000.0;
	tau_ha = 35000.0;

// IP3 Receptor
	//conIP3 = 0.015;			// default ip3 concentration
	//conIP3_inf = 0.015;		// minimum concentration of IP3
	//conIP3_inf = 0.15;		// minimum concentration of IP3
	conIP3_tau = 11000.0;	// usually 1 sec
	//conIP3_cyt = 0.015;		// time dependent change in the cytosolic space
	Pip3_wholeCell = 80.0;	//Maximum Ca permeability of IP3R
	//scfmaxP_IP3R;			// extent of desensitization     given by a scaling factor
	k_1a = 0.64 * 0.001;	// 1/ms   rate constant
	k_1b = 0.04 * 0.001;	// 1/ms   rate constant
	k_2a = 37.4 * 0.001;	// 1/ms   rate constant
	k_2b = 1.4 * 0.001;		// 1/ms   rate constant
	k_3a = 0.11 * 0.001; 	// 1/ms   rate constant
	k_3b = 29.8 * 0.001;	// 1/ms   rate constant
	k_4a = 4 * 0.001;		// 1/ms   rate constant
	k_4b = 0.54 * 0.001;	// 1/ms   rate constant
	l_2a = 1.7 * 0.001;		// 1/ms   rate constant
	l_2b = 0.8 * 0.001;		// 1/ms   rate constant
	l_4a = 1.7 * 0.001;		// 1/ms   rate constant
	l_4b = 2.5 * 0.001;		// 1/ms   rate constant
	l_6a = 4707.0 * 0.001;	// 1/ms   rate constant
	l_6b = 11.4 * 0.001;	// 1/ms   rate constant

	L_1  = 0.12;
	L_3  = 0.025;
	L_5  = 54.7;

// Myoplasmic Ca Ion Concentration Changes, Max. [Ca] buffered in CMDN (mM)  
	KonTnChCa = 2.37;		// 1/(ms)/(mM)
	KoffTnChCa= 0.000032;	// 1/(ms)
	BtotTnCh  = 0.12;		// 0.14 (mM)

	KonCaM  = 34.0;			// 1/(ms)/(mM)
	KoffCaM = 0.238;		// 1/(ms)
	BtotCaM = 0.024;		// mM

	KonSR = 100.0;			// 1/(ms)/(mM)
	KoffSR = 0.06;			// 1/(ms)
	BtotSR = 0.0171;		// 19*0.0009 (mM)

	KonCsqn = 100.0;		// 1/(ms)/(mM)
	KoffCsqn = 65.0;		// 1/(ms)
	Btot_Csqn = 3.0;
	KdCsqnCa = KoffCsqn/KonCsqn;

// Myoplasmic Calcium Buffers
	BtotL_iz = 0.6078;	// 0.0374 * Vol_blk / Vol_iz (mM)
	BtotH_iz = 0.2178;  // 0.0134 * Vol_blk / Vol_iz (mM)

	KoffL_iz = 1.3;		// (ms^-1) sarcolemmal Ca buffer with a low affinity
	KonL_iz  = 100.0;	// (ms^-1mM^-1) sarcolemmal Ca buffer with a low affinity
	KoffH_iz = 0.03;    // (ms^-1)
	KonH_iz  = 100.0;   // (ms^-1mM^-1)
	KdL_iz   = KoffL_iz / KonL_iz;	// dissociation constant of the low affinity buffer in the iz/jnc space
	KdH_iz   = KoffH_iz / KonH_iz;  // dissociation constant of the high affinity buffer in the iz/jnc space

	BtotL_jnc = 1.1095;	// 0.00046*Vol_blk/Vol_jnc (mM) Note, the on- and off-rates are the same as in sl
	BtotH_jnc = 0.398;  // 0.000165*Vol_blk/Vol_jnc (mM)

	GCa_jnciz = 3395.88*Sc;	// 3215.8*0.5*Sc_Vol= 459.4*20*0.7*0.5, 824.13054227789689*ratioCm*ratioVolume
	GCa_izblk = 3507.78*Sc;	// 2076.1*0.8*Sc_Vol= 3724.25607984805*ratioCm*ratioVolume
	
// beta1-adrenergic receptor 
	// total b1-adrenergic receptor, 0.0000132 (mM)
	ARtot = 0.0000132;
	// initial value, total ligand (isoproterenol) (mM)
	//Ligandtot = 0.0001;
	//Ligand = 0.0001;
	//Ligandtot = 0.1;
	Ligand = Ligandtot;
	// total Gs protein, 0.00383 (mM)
	Gstot  = 0.00383;
	// alpha, beta, gamma subunit of Gs (mM)
	Gsfree = Gstot;
	alpha_DS464 = 0.0000011; // ms-1,  forward rate of desensitization
	beta_DS464  = 0.0000022*2.0; // ms-1,  backward  rate of desensitization
	alpha_DS301 = 0.0036; // ms-1,  forward rate of desensitization
	beta_DS301  = 0.0000002232*20.0; // ms-1,  backward  rate of desensitization
	scfalpha_DS464 = 1.0;
	scfalpha_DS301 = 1.0;
	ACtot = 0.0000497; // total adenylate cyclase, 0.0000497 (mM)
	PDE   = 0.000039; // phosphodiesterase, 0.000039 (mM)
	ATPtot = 6.0; // mM
	PKAtot = 0.001;  // total protein kinase A (dimer), 0.001 (mM)
	PKItot = 0.00018;// total protein kinase inhibitor, 0.00018 (mM)
	Cattot = PKAtot*2.0;
	
	alpha = 20.0; // ms-1, an arbitraly value,  A2RC --> Cat + A2R

	// disscociation constant for the instantaneous equilibrium 
    Kd4   = 0.0002;   // mM   PKI + Cat = CPKI
	KdLR  = 0.001;    // mM   R + L = LR,   Saucerman 0.285 micromole
	KdGLR = 0.000062; // mM   LR + G = LRG   Saucerman 0.062 micromole
	KdGR  = 0.033;    // mM    G + R = RG    Saucerman 33 micromole

	KdACG = 0.315;    // AC + GaGTP = AC_GaGTP  'activation of adenylcyclase
	Kd12  = 0.008;    // RC + cAMP = ARC,  ARC + cAMP = A2RC    dissociation constant for the binding of cAMP to R subunit
	Kd3   = 0.009;    // Cat + A2R = A2RC     dissociation constant for the binding of Cat to A2RC

	// action of Cat
	Kf = 0.0625;
	tau_phosphorylation = 40.0*1000.0; // (msec) model adjusted assuming a time constant of 40 sec

// cytosolic substrates for ATP homeostasis
	Mg_cyt = 0.8;

}

