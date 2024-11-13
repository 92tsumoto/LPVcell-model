#include "syspara.h"

void function(double x[],double f[],double t)
{

	// uncheck comp_bAR: 2024/09/19
	comp_aAR_kinetics(x);
	comp_bAR_kinetics(x);
	// check as follows: 2024/09/19
	comp_Ca_buff_blk(x);
	comp_Ca_buff_iz(x);
	comp_Ca_buff_jnc(x);
	comp_buffSR_rl(x);
	comp_bund_dif(x);
	// check as follows: 2024/09/19
	comp_potential(x);
	comp_ina(x);
	comp_ical(x);
	//comp_ical_iz(x);
	//comp_ical_blk(x);
	comp_ikur(x);
	comp_ito(x);
	comp_ik1(x);
	comp_ikr(x);
	comp_inab(x);
	comp_ikb(x);
	comp_icab(x);
	//comp_icab_blk(x);
	//comp_icab_iz(x);
	comp_inak(x);
	comp_INCX(x);
	//comp_INCX_jnc(x);
	//comp_INCX_iz(x);
	//comp_INCX_blk(x);
	comp_ipmca(x);
	//comp_ipmca_iz(x);
	//comp_ipmca_blk(x);
	// check as follows: 2024/09/19
	comp_CaRU(x);	// CaRU
	comp_jtr(x);	// Transfer CaSR
	comp_jup(x);	// SERCA_TC2009_betaadrenergic
	comp_IP3R(x);	// IP3R

	ICaL = ICaLCa_LR + ICaLCa_L0 + ICaLNa_jnc + ICaLK_jnc + ICaLCa_iz + ICaLNa_iz + ICaLK_iz + ICaLCa_blk + ICaLNa_blk + ICaLK_blk;
	ICaL_jnc = ICaLCa_LR + ICaLCa_L0 + ICaLNa_jnc + ICaLK_jnc;
	ICaL_iz = ICaLCa_iz + ICaLNa_iz + ICaLK_iz;
	ICaL_blk = ICaLCa_blk + ICaLNa_blk + ICaLK_blk;
	ICabg_blkiz = ICabg_blk + ICabg_iz;
	totINCX = INCXNa_iz + INCXCa_iz + INCXNa_blk + INCXCa_blk + INCXNa_jnc + INCXCa_jnc;
	totIPMCA = IPMCA_iz + IPMCA_blk;

	Itotal_Ca_jnc = ICaLCa_LR + ICaLCa_L0 - 2.0*INCX_jnc;
	Itotal_Ca_iz  = ICaLCa_iz + ICabg_iz - 2.0*INCX_iz + IPMCA_iz;
	Itotal_Ca_blk = ICaLCa_blk + ICabg_blk - 2.0*INCX_blk + IPMCA_blk;
	Itotal_Ca = Itotal_Ca_jnc + Itotal_Ca_iz + Itotal_Ca_blk;

	Itotal_Na = INaT_Na + (ICaLNa_jnc + ICaLNa_iz + ICaLNa_blk) + INabg + 3.0*INaK + 3.0*(INCX_jnc + INCX_iz + INCX_blk);
	//Itotal_K  = INaT_K + IK1 + iKto + IKur + IKr + IKbg + (ICaLK_jnc + ICaLK_iz + ICaLK_blk) + INaK_K + var.Istim;
	Itotal_K  = INaT_K + IK1 + IKto + IKur + IKr + IKbg + (ICaLK_jnc + ICaLK_iz + ICaLK_blk) - 2.0*INaK;
	
	//printf("ICaLCa_LR=%e,ICaLCa_L0=%e,INCX_jnc=%e,Itotal_Ca=%e\n",ICaLCa_LR,ICaLCa_L0,INCX_jnc,Itotal_Ca);
	//printf("Itotal_Ca_jnc=%e,Itotal_Ca_iz=%e,Itotal_Ca_blk=%e,Itotal_Ca=%e\n",Itotal_Ca_jnc,Itotal_Ca_iz,Itotal_Ca_blk,Itotal_Ca);
	//printf("Inak=%e,Incx_jnc=%e,Incx_iz=%e,Incx_blk=%lf,Inab=%lf,Ina=%lf\n",INaK,INCX_jnc,INCX_iz,INCX_blk,INabg,INaT_Na);
	//printf("Inatk=%e,Ik1=%e,Ikto=%e,Ikur=%lf,Ikr=%lf,Ikbg=%lf\n",INaT_K,IK1,IKto,IKur,IKr,IKbg);
	//printf("Itotal_Na=%e,Itotal_K=%e,Itotal_Ca=%e\n",Itotal_Na,Itotal_K,Itotal_Ca);


// Vm membrane potential
	f[0] = -(Itotal_Ca + Itotal_Na + Itotal_K + var.Istim);
//Fast INa
	f[1] = x[2]*kI2O + fC*C_TM*kC2O - x[1]*(kOC + kOI2); // O_TM for INa
	f[2] = fC*C_TM*kC2I2 + x[1]*kOI2 + x[3]*kIsb - x[2]*(kI2C + kI2O + kIsf); // I2_TM for INa
	f[3] = x[2]*kIsf + C_TM*kIsf - x[3]*2.0*kIsb; // Is_TM for INa
// Ultrarapid IKur
	f[4] = (ass - x[4])/tau_aur;// a_ur
	f[5] = (iss - x[5])/tau_iur; // i_ur
//Transient outward current
	f[6] = (pr_inf - x[6])/tau_pr;
//	f[7] = (s_inf - x[7])/tau_s; // no effect on ito in this model
	f[7] = (ss_inf - x[7])/tau_ss;
//L-type Ca2+ current
	f[8]  = kccco_iz*Ycc_iz + kooco_iz*x[10] - (kcocc_iz + kcooo_iz)*x[8];
	f[9] = koooc_iz*x[10] + kccoc_iz*Ycc_iz - (kocoo_iz + koccc_iz)*x[9];
	f[10] = kcooo_iz*x[8] + kocoo_iz*x[9] - (kooco_iz + koooc_iz)*x[10];
	f[11] = kccco_blk*Ycc_blk + kooco_blk*x[13] - (kcocc_blk + kcooo_blk)*x[11];
	f[12] = koooc_blk*x[13] + kccoc_blk*Ycc_blk - (kocoo_blk + koccc_blk)*x[12];
	f[13] = kcooo_blk*x[11] + kocoo_blk*x[12] - (kooco_blk + koooc_blk)*x[13];
// IKr 
	f[14] = ay1*(1.0 - x[14]) - by1*x[14];
	f[15] = ay2*(1.0 - x[15]) - by2*x[15];
	f[16] = ay3*(1.0 - x[16]) - by3*x[16];
// IK1
	f[17] = bSPM*po_Mggate*(1.0 - x[17]) - aSPM*x[17];
// INCX
	f[18] = E2NCX_jnc*alphaE_jnc + x[19]*beta1_jnc + x[20]*beta2_jnc - x[18]*(betaE_jnc + alpha1_jnc + alpha2_jnc);
	f[19] = x[18]*alpha1_jnc - x[19]*beta1_jnc;
	f[20] = x[18]*alpha2_jnc - x[20]*beta2_jnc;
	f[21] = E2NCX_iz*alphaE_iz + x[22]*beta1_iz + x[23]*beta2_iz - x[21]*(betaE_iz + alpha1_iz + alpha2_iz);
	f[22] = x[21]*alpha1_iz - x[22]*beta1_iz;
	f[23] = x[21]*alpha2_iz - x[23]*beta2_iz;
	f[24] = E2NCX_blk*alphaE_blk + x[25]*beta1_blk + x[26]*beta2_blk - x[24]*(betaE_blk + alpha1_blk + alpha2_blk);
	f[25] = x[24]*alpha1_blk - x[25]*beta1_blk;
	f[26] = x[24]*alpha2_blk - x[26]*beta2_blk;
// IRyR (CaRU)
	f[27] = kYoocYooo*x[28] + kYcooYooo*x[29] + kYocoYooo*x[32] - (kYoooYooc + kYoooYcoo + kYoooYoco)*x[27];
	f[28] = kYcocYooc*x[30] + kYoooYooc*x[27] + kYoccYooc*x[33] - (kYoocYcoc + kYoocYooo + kYoocYocc)*x[28];
	f[29] = kYcocYcoo*x[30] + kYoooYcoo*x[27] + kYccoYcoo*x[31] - (kYcooYcoc + kYcooYooo + kYcooYcco)*x[29];
	f[30] = kYcooYcoc*x[29] + kYoocYcoc*x[28] + kYcccYcoc*(1.0 - (x[27]+x[28]+x[29]+x[30]+x[31]+x[32]+x[33]))  - (kYcocYcoo + kYcocYooc + kYcocYccc)*x[30];
	f[31] = kYcccYcco*(1.0 - (x[27]+x[28]+x[29]+x[30]+x[31]+x[32]+x[33]))  + kYocoYcco*x[32] + kYcooYcco*x[29] - (kYccoYccc + kYccoYoco + kYccoYcoo)*x[31];
	f[32] = kYoccYoco*x[33] + kYccoYoco*x[31] + kYoooYoco*x[27] - (kYocoYocc + kYocoYcco + kYocoYooo)*x[32];
	f[33] = kYcccYocc*(1.0 - (x[27]+x[28]+x[29]+x[30]+x[31]+x[32]+x[33])) + kYocoYocc*x[32] + kYoocYocc*x[28] - (kYoccYccc + kYoccYoco + kYoccYooc)*x[33];
// IP3R model
	f[34] = phi_2b*x[35] + (k_1b + l_2b)*x[36] - (phi_2*conIP3 + phi_1)*x[34];
	f[35] = phi_2*conIP3*x[34] + phi_4b*x[38] + k_3b*(1.0-x[34]-x[35]-x[36]-x[37]-x[38]) - (phi_2b + phi_4 + phi_3)*x[35];
	f[36] = phi_1*x[34] - (k_1b + l_2b)*x[36];
	f[37] = phi_5*x[38] - (k_1b + l_2b)*x[37];
	f[38] = phi_4*x[35] + (k_1b + l_2b)*x[37] - (phi_4b + phi_5)*x[38];
// [Ca]jnc --> Cafree.jnc
	f[39] = Buf_jnc*((-Itotal_Ca_jnc*Cm/(2.0*F) + Jrel_SR - JCa_jnciz + J_ip3R)/Vjnc);
// [Ca]iz --> Cafree.iz
	f[40] = Buf_iz *((-Itotal_Ca_iz*Cm/(2.0*F) + JCa_jnciz - JCa_izblk)/Viz);
// [Ca]blk --> Cafree.blk
	f[41] = Buf_blk*((-Itotal_Ca_blk*Cm/(2.0*F) - JSERCA + JCa_izblk) / Vblk);
// [Ca]SR_uptake--> Ca_SRup
	f[42] = (JSERCA - Jtrans_SR) / Vsr_up;
// [Ca]SR_release--> Cafree_SRrl
	f[43] = Buf_SRrl*((Jtrans_SR - Jrel_SR - J_ip3R) / Vsr_rl);
// [Na]i--> Nai
	f[44] = -Itotal_Na*Cm/(Vcyt*F);
// [K]i --> Ki
	f[45] = -(Itotal_K + var.Istim)*Cm/(Vcyt*F);
// bAR model
	// ARS464
	//f[48] = alpha_DS464*scfalpha_DS464*(LAR + LARGs) - beta_DS464*x[48]; 
	f[46] = alpha_DS464*(LAR + LARGs) - beta_DS464*x[46]; 
	// ARS301
	//f[49] = alpha_DS301*scfalpha_DS301*x[54]*(ARtot - x[48] - x[49]) - beta_DS301*x[49];
	f[47] = alpha_DS301*x[52]*(ARtot - x[46] - x[47]) - beta_DS301*x[47];
	// GsaGTP
	f[48] = B2*(0.016*(ARGs + LARGs) - 0.001*(x[48] + ACtot*x[48]/(x[48] + KdACG)));
	// Gsbg
	f[49] = 0.016*(ARGs + LARGs) - 1200.0*x[49]*x[50];
	// GsaGDP
	f[50] = 0.001*(x[48] + ACtot*x[48]/(x[48] + KdACG)) - 1200.0*x[49]*x[50];
	// cAMP
	f[51] = B1*(0.0001307*AC*ATPtot/(1.03 + ATPtot) + 3.4*AC_GaGTP*ATPtot/(0.315 + ATPtot) - 0.005*PDE*x[51]/(0.0013 + x[51]));
	// Cat (cataliticC)
	f[52] = B0*(beta_PKA*(2.0*PKAtot) - (alpha_PKA + beta_PKA)*(x[52] + PKItot*x[52]/(x[52] + Kd4)));
	// y_PKA
	//f[53] = kFcat*(1.0 - x[53]) - kBpp*x[53];
	f[53] = kFcat*(1.0 - x[53]) - kBpp*x[53];
/*
	//f[0]=0.0;
	//f[1]=0.0;
	//f[2]=0.0;
	//f[3]=0.0;
	f[4]=0.0;
	f[5]=0.0;
	f[6]=0.0;
	f[7]=0.0;
	f[8]=0.0;
	f[9]=0.0;
	f[10]=0.0;
	f[11]=0.0;
	f[12]=0.0;
	f[13]=0.0;
	f[14]=0.0;
	f[15]=0.0;
	f[16]=0.0;
	f[17]=0.0;
	f[18]=0.0;
	f[19]=0.0;
	f[20]=0.0;
	f[21]=0.0;
	f[22]=0.0;
	f[23]=0.0;
	f[24]=0.0;
	f[25]=0.0;
	f[26]=0.0;
	f[27]=0.0;
	f[28]=0.0;
	f[29]=0.0;
	f[30]=0.0;
	f[31]=0.0;
	f[32]=0.0;
	f[33]=0.0;
	f[34]=0.0;
	f[35]=0.0;
	f[36]=0.0;
	f[37]=0.0;
	f[38]=0.0;
	f[39]=0.0;
	f[40]=0.0;
	f[41]=0.0;
	f[42]=0.0;
	f[43]=0.0;
	f[44]=0.0;
	f[45]=0.0;
	f[46]=0.0;
	f[47]=0.0;
	f[48]=0.0;
	f[49]=0.0;
	f[50]=0.0;
	f[51]=0.0;
	f[52]=0.0;
	f[53]=0.0;
	*/
}
