//#ifndef __SYSPARA_H_INCLUDE 
//#define __SYSPARA_H_INCLUDE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mkl.h"
#include "/home/tsumoto/lib/icx2024/xhplot.h"

#define NN 54
#define BUF 100
#define NUM 20

#define R 8.3143	// C*mV/mmol/K
#define F 96.4867	// C/mmol
#define T 310.0		// K

#define dvm 5
#define Emax 2000
#define Emin -2000
#define VNMAX (Emax-Emin)*dvm+1

	// Cell geometry
	//double length;	// cell length
	double Sc;          // scaling factor inducing from HuVEC model
	double Vcell;       // Cell Volume (pL): 120x37.62x8.4xscaling-factor
    double Vblk;		// Bulk(myoplasm) volume (pL)
	double Viz;			// IZ volume (pL)
    double Vjnc;		// JSR volume (pL)
	double Vcyt; 		// cytosolic space (pF)
	double Vsr;			// whole SR volume (pL)
	double Vsr_rl;		// fractional volume of the junctional SR (pL)
    double Vsr_up;		// fractional volume of the network SR (pL)
	
	double Cm;			// 113.5514 (pF):Capacitive membrane area (pF)
	// only for the mitochondria model
	double Rcm_mit;		// volume ratio of cytosol/mitochondria (Vol_cyt / Vol_mit):
	double Vmit;		// Mitochondria volume (pL)
	double Cm_mit;		// Mitochondria surface membrane capacitance (pF)
	double ampRatio;	// Vcyt * vRatio_mit <-- ? From Umehara et al., IJMS, 2019

	// Fraction of Currents jnc, iz, blk
	double Frc_iz;		// fractional distribution at the iz myoplasm for ion transporters except LCC and NCX
    double Frc_blk;		// fractional distribution at the blk myoplasm for ion transporters except LCC and NCX
	double Frc_CaL_jnc;	// fractional distribution of ICaL at the junctional myoplasm
	double Frc_CaL_iz;	// fractional distribution of ICaL at the iz myoplasm
	double Frc_CaL_blk;	// fractional distribution of ICaL at the blk myoplasm
	double Frc_NCX_jnc;	// fractional distribution of NCX at the junctional myoplasm
	double Frc_NCX_iz;	// fractional distribution of NCX at the iz myoplasm
	double Frc_NCX_blk;	// fractional distribution of NCX at the blk myoplasm

	// Ion Valences
	double zna,zk,zca; 

	// An invariant constant
	double RTonF,RTon2F;

	// Extracellular Ion concentrations
	double Nao,Cao,Ko; //nai;

	// Other Ion concentrations
	double Mg_cyt;
	double Pifree_cyt,AMP_cyt,ADPfr_cyt,MgADP_cyt;
	double ATPfr_cyt,MgATP_cyt,Cr_cyt,H_cyt;

	// Reversal potential
	double Ena,Ek,Eks;
	double ECa_jnc,ECa_iz,ECa_blk;
	double ECa_LR,ECa_L0;
	double GHK_Ca_jnc,GHK_Ca_iz,GHK_Ca_blk;
	double GHK_Ca_LR,GHK_Ca_L0;
	double GHK_Na,GHK_K,GHK_Ca;

	// Ion currents 
	double Itotal_Ca,Itotal_Na,Itotal_K;
	double INaT_Na,INaT_K, INaT;
	double IK_TM_K_cyt,  Ik_ss_total;
	double Ica_i_total,Ica_ss_total;
	double Itotal_Ca_jnc,Itotal_Ca_iz,Itotal_Ca_blk;
	double Itotal_Na_cyt,Itotal_K_cyt;
	// Base Currnt Stimulus
	double Istim;

struct varstruct {

    int datas;
    int line_wid[NUM];
	
	int n;
    double Istim;
	double dIstim;

	// Base Currnt Stimulus
	double Istim_base;
	double omega,ndis;

	// test variable
	double dt;
	// Sttimulus parameters
	double BCL;  // Base cycle length = stimulus period
	int beat; // Number of stimulus

    int m;
    int l;

    //double x0[NUM][NN];
    double **x0;
    double tsign[NUM];
    double tend[NUM];

    int pflag;
    int write, graph;
    int write0;
    int half;
    int debug;
	int pswitch;

} var;

// time-dependent active fraction of target protein phosppholylated
	double afLCC,baseCaL,deltaCaL; // LCC phosppholylated by PKA + base (activated by some other PKs)
	double afSERCA,baseSC,deltaSC; // phospholamban of SERCA phospholylated by PKA + base (activated by some other PKs)
	double afNaK,baseNaK,deltaNaK; // phospholemman of NaK phospholylated by PKA + base (activated by some other PKs)
	double ATPfactor;

// Fast and Late sodium currnets
	double *TkC2O,*TkOC,*TkOI2,*TkIsb,*TkIsf,*TfC; // Himeno et al., PVC model
	double kC2O,kOC,kOI2,kIsb,kIsf,C_TM,fC,kC2I2,kI2C; // Himeno et al., PVC model
	double kI2O,kI1I2,kI1O,fL;
	double PNa, rPK_Ina,rate_INa;

// Transient Outward Current (Ito)
	double IKto,GKto,rate_Ito,fast_ito,slow_ito;
	double *Tpr_inf,*Ttau_pr,*Ts_inf,*Ttau_s,*Tss_inf,*Ttau_ss;
	double pr_inf,tau_pr,s_inf,tau_s,ss_inf,tau_ss;

// L-type Calcium channel current (IcaL)
	double ICaL;
	double ICaLCa_LR,ICaLCa_L0;
	double ICaL_jnc,ICaL_iz,ICaL_blk;
	double ICaLCa_iz,ICaLCa_blk;
	double ICaLNa_jnc,ICaLNa_iz,ICaLNa_blk;
	double ICaLK_jnc,ICaLK_iz,ICaLK_blk;
	double KL,TL,gD_free,JLCC;
	double Ycc_iz, kcocc_iz, kcooo_iz, kccco_iz, kooco_iz, koooc_iz, kocoo_iz, koccc_iz, kccoc_iz;
	double Ycc_blk,kcocc_blk,kcooo_blk,kccco_blk,kooco_blk,koooc_blk,kocoo_blk,koccc_blk,kccoc_blk;
	double VshiftLCCact,VshiftLCCinact;
	double *TfVmACT,*TexpF;
	double *Tkooco_iz,*Tkccco_iz;
	double *Tkooco_blk,*Tkccco_blk;
	double scf_LCC,tPCaL,oPCaL,ATPCaL;

// Rapid activating potassium current (Ikr)
	double IKr_K_cyt,IKr,rate_Kr,gKr,qKr;
	double *Tay1,*Tby1,*Tay2,*Tby2,*Tay3,*Tby3;
	double ay1,by1,ay2,by2,ay3,by3;

// Ultrarapidactivating potassium current (Ikur)
	double IKur,GKur;
	double *Tass,*Tiss,*Ttau_aur,*Ttau_iur;
	double ass,iss,tau_aur,tau_iur;
	double rate_Kur;

// Cloride current from Okamoto et al., JMCC, 2012,2014
	double *Ta1_Clh,*Ta2_Clh,*Tb1_Clh,*Tb2_Clh;
	double a1_Clh,a2_Clh,b1_Clh,b2_Clh;
	double IClh,rate_Clh,GClh;

// Inward rectifier potassium current (Ik1)
	double IK1_K_cyt,IK1,frac_mode_K1,GK1;
	double aSPM,bSPM,po_Mggate;
	double SPM,qK1,fK1;

// Sodium-Calcium exchanger (NCX)
	double rate_NCX_jnc,rate_NCX_iz,rate_NCX_blk;
	double maxINCX;
	double totINCX;
	double INCXCa_jnc,INCXCa_iz,INCXCa_blk;
	double INCXNa_jnc,INCXNa_iz,INCXNa_blk;
	double INCX_jnc,INCX_iz,INCX_blk;
	double *Tk1_NCX,*Tk2_NCX;
	double a1on,a1off,a2on,a2off;
	double b1on,b1off,b2on,b2off;
	double KmNao,KmNai,KmCao,KmCai;
	double Kmact,k3,k4,E2Na,E2Ca;
//	NCX_jnc
	double alpha1_jnc,alpha2_jnc,alphaE_jnc;
	double beta1_jnc,beta2_jnc,betaE_jnc;
	double E2NCX_jnc;
//	NCX_iz
	double alpha1_iz,alpha2_iz,alphaE_iz;
	double beta1_iz,beta2_iz,betaE_iz;
	double E2NCX_iz;
//	NCX_blk
	double alpha1_blk,alpha2_blk,alphaE_blk;
	double beta1_blk,beta2_blk,betaE_blk;
	double E2NCX_blk;

// Na-K Pump
	double INaK,INaK_Na,INaK_K;
	double k1p,k1m,k2p,k2m,k3p,k3m,k4p,k4m;
	double KdKe0,KdKi0,KdNae0,KdNai0,KdMgATP;
	double stoiNa,stoiK;

	double dV_Nae,dV_Nai,dV_Ke,dV_Ki;
	//double sfKdNaK;    // modulation by the phosphlemman

	double maxINaK,rate_NaK;

// Sarcolemmal Ca Pump
	double totIPMCA;
	double IPMCA_iz,IPMCA_blk;
	double maxIPMCA,KmPMCA;

// Na Background Current
	double INabg,gnab;

// K Background Current
	double IKbg,gkb;
	double scfIKb,tauIKbg; 

// Ca Background Current
	double ICabg_iz,ICabg_blk,ICabg_blkiz;
	double ICabg,gcab;

// SR calcium release flux, via RyR (Jrel)
	double Jrel_SR,PRyR;
	double kYoooYooc,kYoooYcoo,kYoooYoco;
	double kYoocYooo,kYoocYcoc,kYoocYocc;
	double kYocoYooo,kYocoYcco,kYocoYocc;
	double kYoccYooc,kYoccYoco,kYoccYccc;
	double kYcooYooo,kYcooYcoc,kYcooYcco;
	double kYcocYooc,kYcocYcoo,kYcocYccc;
	double kYccoYcoo,kYccoYccc,kYccoYoco;
	double kYcccYocc,kYcccYcco,kYcccYcoc;
	double Yccc;
	double JRyR,JLCC,gD_nd;
	double Ca_Noise,Ca_ndLR,Ca_ndL0;
	double *Tkvoc,*Tkcaco;
	double Jleak_SR;	// not use?

// Calcium uptake via SERCA pump
	double JSERCA;
	double maxISERCA;
	double Kd_H1,Kd_Hi,Kd_Hsr,Kd_H_release;
	double Hi,Hsr;

// Translocation of Ca Ions from NSR to JSR
	double Jtrans_SR,Ptrans;
	double Cafree_SRrl,Buf_SRrl;

//	IP3 Receptor 
	double J_ip3R;
	double conIP3,conIP3_inf,conIP3_tau,conIP3_cyt,Pip3_wholeCell,scfmaxP_IP3;
	double k_1a,k_1b,k_2a,k_2b,k_3a,k_3b,k_4a,k_4b;
	double l_2a,l_2b,l_4a,l_4b,l_6a,l_6b;
	double L_1,L_3,L_5;
	//double uMCa;	// Ca concentration microM
	double phi_1,phi_2,phi_2b;
	double phi_3,phi_4,phi_4b;
	double phi_5;
	double S_ip3;
// aAR
	double m_inf,h_inf;
	double tau_mi,tau_ha;

// bAR
	double kFcat,kBpp;
	// total b1-adrenergic receptor, 0.0000132 (mM)
	double ARtot;
	// initial value, total ligand (isoproterenol) (mM)
	double Ligandtot;
	double Ligand; // assumed to be = Ligandtot in the present simulation   (mM))
	// free b1-adrenoceptor (mM)
	double ARfree;
	double ARact; //active b1-adrenergic receptor (mM)
	// b1-adrenergic receptor bound with ligand (mM)
	double LAR;
	// b1-adrenergic receptor bound with ligand and Gs (mM)
	double LARGs;
	// b1-adrenergic receptor bound with Gs (mM)
	double ARGs;
	// total Gs protein, 0.00383 (mM)
	double Gstot;
	// alpha, beta, gamma subunit of Gs (mM)
	double Gsfree;
	double alpha_DS464,beta_DS464,alpha_DS301,beta_DS301;
	double scfalpha_DS464,scfalpha_DS301;
	// total adenylate cyclase, 0.0000497 (mM)
	double ACtot;
	// temporary buffering coefficent
	double AAC;
	// adenylate cyclase (mM)
	double AC;
	// adenylate cyclase-bound GsaGTP (mM)
	double AC_GaGTP;
	// phosphodiesterase, 0.000039 (mM)
	double PDE;
	double ATPtot;
	double B2;
	// free cAMP
	double B1,B10,B11,B12,B13;
	double PKAtot;// total protein kinase A (dimer), 0.001 (mM)
	// regulatory subunit-catalytic subunit complex.  (mM)
	double RC;
	// complex bound with one cAMP (mM)
	double ARC_bAR;
	// complex bound with 2 cAMP (mM)
	double A2RC;
	// regulatory subunit bound with 2 cAMP (mM)
	double A2R;
	//double cataliticC; // catalytic subunit of protein kinase A
	double B0; // catalytic subunit of protein kinase A
	double C_PKI; // PKA - PKI complex (mM)

	double PKI; // protein kinase inhibitor (mM)
	double PKItot; // total protein kinase inhibitor, 0.00018 (mM)
	double alpha_PKA;// ms-1, Cat123 --> Cat45
	double beta_PKA; // ms-1, Cat123 <-- Cat45
	double Cat123;   // mM
	double Cattot;
	double Cat45;
	double ACAT;

	//double GsaGTPtot; // total GTP-bound alpha subunit of Gs (mM)
	double alpha,beta;
	
	//disscociation constant for the instantaneous equilibrium
	double Kd4,KdLR,KdGLR,KdGR,KdACG,Kd12,Kd3;

	// action of Cat
	double Kf,kFcat; // Vmax, Kf*[catalyticC]
	double tau_phosphorylation; // (msec): model adjusted assuming a time constant of 40 sec



// Myoplasmic Calcium Buffers
	double BtotL_iz,BtotH_iz;
	double KoffL_iz,KonL_iz,KoffH_iz,KonH_iz;
	double KdL_iz,KdH_iz;
	double Buf_iz;

	double BtotL_jnc,BtotH_jnc;
	double Buf_jnc;
	double Buf_blk;
	double BtotL_blk,BtotH_blk;
	double KonTnChCa,KoffTnChCa,BtotTnCh;
	double KonCaM,KoffCaM,BtotCaM;
	double KonSR,KoffSR,BtotSR;

	double KonCsqn,KoffCsqn,Btot_Csqn,KdCsqnCa;

	double JCa_jnciz,GCa_jnciz;
	double JCa_izblk,GCa_izblk;

	double KdCaM;
	double dCafree_blkdt;
	double freezeSR;

//int main(int argc,char **argv);
void val_consts(FILE *);
void make_ExpTable(void);

//void eular(int n,double h,double x[],double t);
void runge(int n,double h,double x[],double t);
void eular(int n,double h,double x[],double t);
void function(double x[],double f[],double t);
void input_para(FILE *);
void initial_mem(void);
void closed_mem(void);

void data_out(FILE *,double t,double x[]);
void current(FILE *,FILE *,FILE *,FILE *,FILE *,FILE *,FILE *,double t,double x[]);

void eventloop(FILE *, int *mode, int *P, double m[]);
void orbit(int *mode, double m[], double x2);
void draw_p(int *mode, int P, double x[], double x2);
void mouse(int *mode, double x[], double x2);

void comp_aAR_kinetics (double x[]);
void comp_bAR_kinetics (double x[]);
void comp_Ca_buff_blk (double x[]);
void comp_Ca_buff_iz (double x[]);
void comp_Ca_buff_jnc (double x[]);
void comp_buffSR_rl (double x[]);
void comp_bund_dif (double x[]);
void comp_potential(double x[]);
void comp_ina(double x[]);
void comp_ical(double x[]);
//void comp_ical_iz(double x[]);
void comp_ical_iz_gating(double x[]);
void comp_ical_iz_amp(double x[]);
void comp_ical_blk(double x[]);
void comp_ical_blk_gating(double x[]);
void comp_ical_blk_amp(double x[]);
void comp_ikur(double x[]);
void comp_ito(double x[]);
void comp_ik1(double x[]);
void comp_ikr(double x[]);
//void comp_iclh(double x[]);
void comp_inab(double x[]);
void comp_ikb(double x[]);
void comp_icab(double x[]);
//void comp_icab_blk(double x[]);
//void comp_icab_iz(double x[]);
void comp_inak(double x[]);
void comp_INCX(double x[]);
//void comp_INCX_jnc(double x[]);
//void comp_INCX_iz(double x[]);
//void comp_INCX_blk(double x[]);
void comp_ipmca(double x[]);
//void comp_ipmca_iz(double x[]);
//void comp_ipmca_blk(double x[]);
void comp_jtr (double x[]);
void comp_CaRU (double x[]);
void comp_jup(double x[]);
void comp_IP3R (double x[]);

void current_ikr(FILE *, double t, double x[]);
void current_ina(FILE *, double t, double x[]);
void current_ical(FILE *, double t, double x[]);
void current_inaca(FILE *, double t, double x[]);
void current_inak(FILE *, double t, double x[]);
void current_ikur(FILE *, double t, double x[]);
void current_ikto(FILE *, double t, double x[]);
//void comp_iconcent (double x[]);
//void comp_iconcent2 (double x[]);
//void conc_nsr(double x[]);
//void conc_jsr(double x[]);
//void conc_itr (double x[]);
//void conc_cai (double x[]);
//void conc_cleft (double x[]);
void *malloc2d(size_t size, int row, int col);
void *malloc3d(size_t size, int i, int j, int k);

