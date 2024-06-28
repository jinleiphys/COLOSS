// Calculation of a Wronskian testing the quality of calculated Coulomb wave functions at a given z.
// -------------------------------------------------------------------------------------------------
//
// One uses the Wronskian F[l].H[omega,l+1] - F[l+1].H[omega,l] = C(l,eta)/C(l+1,eta)/(2l+3) and its derivative, which is zero, with standard normalization.
// One calculates |H[omega,l+1]/H[omega,l] - F[l+1]/F[l] - C(l,eta)/C(l+1,eta)/(2l+3)/(F[l].H[omega,l])|oo
// and |(H[omega,l+1]/H[omega,l].(F'[l]/F[l]) + H'[omega,l+1]/H[omega,l] - F'[l+1]/F[l] - (H[omega,l+1]/H[omega,l]).(F[l+1]/F[l])|oo
// omega is chosen so H[omega] is minimal if F is not, i.e. |H[omega,l]| <= |H[-omega,l]|.
// As recurrence relations are not used in the class Coulomb_wave_function, this is a meaningful test.
// The other Wronskian F'.H[omega] - F.H[omega'] = 1 does not test anything,
// as it is used with continued fractions to evaluate Coulomb wave functions.
// For non-standard normalization, one has to use instead F[l].H[omega,l+1] - F[l+1].H[omega,l].[(C(l+1,eta)/C(l,eta))^2] = 1/(2l+3)
// 
// Variables
// ---------
// z : variable of the wave functions
// cwf_l,cwf_lp1,l,eta : classes Coulomb_wave_function of respective parameters l,eta and l+1,eta for angular momentum and Sommerfeld parameter respectively. 
//                       Normalization can be standard or not, but it has to be the same for both classes.
// is_it_normalized : true if standard normalization is used, false if not.
// Fl,dFl,Flp1,dFlp1 : regular wave functions and derivative calculated at z of parameters l,eta and l+1,eta respectively.
// Hp_l,Hm_l,dHp_l,dHm_l : irregular wave functions H+, H- and derivatives calculated at z of parameters l and eta.
// omega: +1 or -1, so F and H[omega] are numerically linearly independent. One then has |H[omega](z)| <= |H[-omega](z)|.
// H_omega_l,dH_omega_l,H_omega_lp1,dH_omega_lp1 : irregular wave functions H+(omega = 1) or H-(omega = -1) and derivative 
//                                                 calculated at z of respective parameters l,eta and l+1,eta.
// Cl_eta_ratio, Cl_eta_ratio_inv_square : C(l,eta)/C(l+1,eta), [C(l+1,eta)/C(l,eta)]^2.
// F_ratio, H_omega_ratio,Fl_log_der,H_omega_l_log_der : F[l+1]/F[l], H[omega,l+1]/H[omega,l], F'[l]/F[l], H'[omega,l]/H[omega,l].
// inf_norm_W_diff : test associated to  F[l], H[omega,l], F[l+1] and H[omega,l+1] (see above).
// inf_norm_W_der_diff :  test associated to  F[l], H[omega,l], F[l+1], H[omega,l+1] and their derivative (see above).

double Wronskian_test (const complex<double> &z,
		       class Coulomb_wave_functions &cwf_l,
		       class Coulomb_wave_functions &cwf_lp1)
{
  const bool is_it_normalized = cwf_l.is_it_normalized;
  const complex<double> l = cwf_l.l, eta = cwf_l.eta;

  if ((cwf_lp1.is_it_normalized != is_it_normalized) || (cwf_lp1.l != l + 1) || (cwf_lp1.eta != eta)) 
    cout<<"Parameters of cwf_l and cwf_lp1 are not correct."<<endl, exit (1);

  complex<double> Fl,dFl,Flp1,dFlp1,Hp_l,dHp_l,Hm_l,dHm_l;

  cwf_l.F_dF (z,Fl,dFl);
  cwf_lp1.F_dF (z,Flp1,dFlp1);
  
  cwf_l.H_dH (1,z,Hp_l,dHp_l);
  cwf_l.H_dH (-1,z,Hm_l,dHm_l);

  const int omega = (abs (Hp_l) <= abs (Hm_l)) ? (1) : (-1);
  const complex<double> H_omega_l = (omega == 1) ? (Hp_l) : (Hm_l), dH_omega_l = (omega == 1) ? (dHp_l) : (dHm_l);

  complex<double> H_omega_lp1,dH_omega_lp1;
  cwf_lp1.H_dH (omega,z,H_omega_lp1,dH_omega_lp1);

  const complex<double> Cl_eta_ratio = exp (log_Cl_eta_calc (l,eta) - log_Cl_eta_calc (l+1,eta)), F_ratio =  Flp1/Fl, H_omega_ratio = H_omega_lp1/H_omega_l;

  if (is_it_normalized)
  {
    const double inf_norm_W_diff = inf_norm (H_omega_ratio - F_ratio - Cl_eta_ratio/(H_omega_l*Fl*(2*l + 3)));
    const double inf_norm_W_der_diff = inf_norm ((dFl/Fl)*H_omega_ratio + dH_omega_lp1/H_omega_l - dFlp1/Fl - (dH_omega_l/H_omega_l)*F_ratio);

    return max (inf_norm_W_diff,inf_norm_W_der_diff);
  }
  else
  {
    const complex<double> Cl_eta_ratio_inv_square = 1.0/(Cl_eta_ratio*Cl_eta_ratio), Fl_log_der = dFl/Fl, H_omega_l_log_der = dH_omega_l/H_omega_l;
    const double inf_norm_W_diff = inf_norm (H_omega_ratio - F_ratio*Cl_eta_ratio_inv_square - 1.0/(H_omega_l*Fl*(2*l + 3)));
    const double inf_norm_W_der_diff = inf_norm (Fl_log_der*H_omega_ratio+dH_omega_lp1/H_omega_l - (dFlp1/Fl+H_omega_l_log_der*F_ratio)*Cl_eta_ratio_inv_square);

    return max (inf_norm_W_diff,inf_norm_W_der_diff);
  }
}







// Calculation of tables of integer l's spaced of Coulomb wave functions with three-term recurrence relations
// ----------------------------------------------------------------------------------------------------------
// 
// Coulomb wave functions obey three-term recurrence relations:
// 
// For f being F,G,H+ or H-, with the standard normalization, one has:
//
// f[l] = (S[l]/R[l]).f[l-1] - (1/R[l]).f'[l-1], f'[l] = R[l].f[l-1] - S[l].f[l]  (forward recurrence) 
// f[l] = (1/R[l+1]).f'[l+1] + (S[l+1]/R[l+1]).f[l+1], f'[l] = S[l+1].f[l] - R[l+1].f[l+1]  (backward recurrence) 
// 
// With the non-standard normalization, one has, for F and f being G,H+ or H- only:
// F[l] = (2l+1).(S[l]/(R[l]^2)).F[l-1] - ((2l+1)/(R[l]^2)).F'[l-1], F'[l] = (2l+1).F[l-1] - S[l].F[l]  (forward recurrence) 
// F[l] = (1/(2l+3]).F'[l+1] + (S[l+1]/2(l+3)).F[l+1], F'[l] = S[l+1].F[l] - ((R[l+1]^2)/(2l+3)).F[l+1]  (backward recurrence)
// 
// f[l] = (S[l+1]/(2l+1]).f[l-1] - (1/(2l+1)).f'[l-1], f'[l] = ((R[l]^2)/(2l+1)).f[l-1] - S[l].f[l]  (forward recurrence) 
// f[l] = ((2l+3)/(R[l+1]^2)).f'[l+1] + (2l+3).(S[l+1]/(R[l+1]^2)).f[l+1], f'[l] = S[l+1].f[l] - (2l+3).f[l+1]  (backward recurrence)  
// 
// where R[l] = (2l+1).C(l,eta)/C(l-1,eta), S[l] = l/z + eta/l .
// l varies from l_deb to l_end = l_deb + Nl - 1.
//
// Forward recurrence is stable for |f| increasing and backward recurrence is stable for |f| decreasing.
// Asymptotically, |F| -> 0 and |H[omega]| -> +oo with Re[l] (the choice of omega is given afterwards).
// Then, forward recurrence is used for H[omega] and backward recurrence for F.
// There may be, however, a turning point l_turn at low |l|, for which |F| increases and |H[omega]| decreases.
// In this case, one has to reverse recurrences for F and H[omega].
// Omega is chosen so |H[omega,l_turn]| <= |H[-omega,l_turn]|.
// With this omega, H[omega] will be minimal for low |l| if F is not minimal in this region.
//
// In recurrence routines, all functions F, G, H+ and H- and derivatives are calculated, 
// with the standard or non-standard normalization, for a single z or a table z_tab of Nz z values.
// f/f' is one-dimensional table for a single z: f(l,eta,z) = f[il], where l = l_deb + il.
// f/f' is a two-dimensional table for a table of z: f(l,eta,z) = f[iz][il], where l = l_deb + il and z = z_tab[iz].
// There also, it is better to have |F[iz][l_end]| with iz increasing.
// Scaling of H+/H- is not included here, as it is useless due to the necessary knowledge of F (i.e., |F| would be out of range if it were.).
// Functions ending with "helper" are made to be called from other functions and do not have to be used by the user.



// Routine calculating all F,F' tables with recurrence relations and the angular turning point if any.
// ---------------------------------------------------------------------------------------------------
// 
// F,F' are calculated for all z's and l's with recurrence relations from the direct evaluation of F,F' at l=l_end..
// l_turn, functions of z, is determined from the condition |F[l_turn]| < |F[l_turn+1]| for standard normalization,
// |F[l_turn]| < |F[l_turn+1].R[l+1]/(2l+3)| for non-standard normalization.
// If there is a turning point, F,F' are directly calculated at l=l_deb and forward recurrence is used.
//
// Variables
// ---------
// l_deb,l_end : first and last considered angular momenta.
// Nl,l : number of angular momenta, angular momentum l = l_deb+i, i in [0:Nl-1]. l_end is l_deb + Nl - 1.
// eta : Sommerfeld parameter
// is_it_normalized : true if one uses standard normalization, false if not.
// Nz : number of z variables
// z_tab : table of Nz variables
// R_tab : Table of Rl values : Rl = (2l+1).C(l,eta)/C(l-1,eta), l in [l_deb+1,l_end].
// cwf_l_deb,cwf_l_end : classes Coulomb_wave_functions of respective parameters l_deb, eta, is_it_normalized and l_end, eta, is_it_normalized.
// F_tab,dF_tab : table of values of F(l,eta,z) and F'(l,eta,z). F(l,eta,z) = F_tab[iz][il], F'(l,eta,z) = dF_tab[iz][il], with z = z_tab[iz] and l = l_deb+il.
// il_turn_tab : l_turn = l_deb + il_turn_tab[iz] if there is a turning l point for a given z (see above).
// il_turn : reference on il_turn_tab[iz].
// z,Fz_tab,dFz_tab : z_tab[iz], F_tab[iz] and dF_tab[iz] for iz in [0:Nz-1].
// Fl,dFl,Flp1,dFlp1,Flm1,dFlm1 : F_tab[iz][il], dF_tab[iz][il], 
//                                F_tab[iz][il+1], dF_tab[iz][il+1],
//                                F_tab[iz][il-1], dF_tab[iz][il-1], respectively.
// Sl,Rl : Sl = S[l] = l/z + eta/l and Rl = R[l] = R_tab[il], for l = l_deb+il and il in [1:Nl-1]. 
// Slp1,Rlp1 : Slp1 = S[l+1] = (l+1)/z + eta/(l+1) and Rlp+1 = R[l+1] = R_tab[il+1], for l = l_deb+il and il in [0:Nl-2].

void F_dF_l_tables_rec_rel_helper (const int Nl,const int Nz,const complex<double> z_tab[],const complex<double> R_tab[],
				   class Coulomb_wave_functions &cwf_l_deb,
				   class Coulomb_wave_functions &cwf_l_end,
				   complex<double> *const *const F_tab,complex<double> *const *const dF_tab,int il_turn_tab[])
{
  const bool is_it_normalized = cwf_l_deb.is_it_normalized;
  const complex<double> l_deb = cwf_l_deb.l,eta = cwf_l_deb.eta;

  for (int iz = 0 ; iz < Nz ; iz++)
  { 
    int &il_turn = il_turn_tab[iz];
    il_turn = -1;

    const complex<double> z = z_tab[iz];
    complex<double> *const Fz_tab = F_tab[iz], *const dFz_tab = dF_tab[iz];

    cwf_l_end.F_dF (z,Fz_tab[Nl-1],dFz_tab[Nl-1]);

    for (int il = Nl-2 ; il >= 0 ; il--)
    {
      const complex<double> l = l_deb+il,Slp1 = (l+1)/z + eta/(l+1),Rlp1 = R_tab[il+1],Flp1 = Fz_tab[il+1],dFlp1 = dFz_tab[il+1];
      complex<double> &Fl = Fz_tab[il],&dFl = dFz_tab[il];

      Fl = (is_it_normalized) ? ((dFlp1 + Slp1*Flp1)/Rlp1) : ((dFlp1 + Slp1*Flp1)/(2*l+3));
      dFl = (is_it_normalized) ? (Slp1*Fl - Rlp1*Flp1) : (Slp1*Fl - (Rlp1*Rlp1/(2*l + 3))*Flp1);

      if (is_it_normalized && (abs (Fl) < abs (Flp1))) 
	il_turn = il, il = -1;
      else if (!is_it_normalized && (abs (Fl) < abs (Flp1*Rlp1/(2*l+3)))) 
	il_turn = il, il = -1;
    }
    
    if (il_turn >= 0)
    {
      cwf_l_deb.F_dF (z,Fz_tab[0],dFz_tab[0]);
      for (int il = 1 ; il <= il_turn ; il++)
      {
	const complex<double> l = l_deb+il,Sl = l/z + eta/l,Rl = R_tab[il],Flm1 = Fz_tab[il-1],dFlm1 = dFz_tab[il-1];
	complex<double> &Fl = Fz_tab[il],&dFl = dFz_tab[il];

	Fl = (is_it_normalized) ? ((Sl*Flm1 - dFlm1)/Rl) : ((Sl*Flm1 - dFlm1)*(2*l+1)/(Rl*Rl));
	dFl = (is_it_normalized) ? (Rl*Flm1 - Sl*Fl) : ((2*l+1)*Flm1 - Sl*Fl);
      }
    }
  }
}







// Routine calculating G,H+,H- and derivative tables for a given z with recurrence relations.
// ------------------------------------------------------------------------------------------
// 
// z is fixed for this routine.
// H+,H+' and H-,H-' are calculated first directly at l=l_turn to determine omega so H[omega,l_turn] is minimal if F is not.
// Then H[omega],H'[omega] are calculated for all l with recurrence relations with forward and backward recurrence from l_turn.
// G,G' and H[-omega],H'-[omega] come forward with standard relations.
//
// Variables
// ---------
// l_deb,l_end : first and last considered angular momenta.
// Nl,l : number of angular momenta, angular momentum l = l_deb+i, i in [0:Nl-1]. l_end is l_deb + Nl - 1.
// is_it_normalized : true if one uses standard normalization, false if not.
// R_tab : table of R[l] values : R[l] = (2l+1).C(l,eta)/C(l-1,eta), l in [l_deb+1,l_end].
// cwf_l_tab : table of pointers to Coulomb_wave_functions classes of respective parameters l=l_deb+il, eta, is_it_normalized for il in [0:Nl-1].
// z,Fz_tab,dFz_tab : parameter of the Coulomb wave functions, tables of F,F' values for a fixed z and varying l.
// Hp_z_tab,dHp_z_tab : table of values of H+(l,eta,z) and H+'(l,eta,z).
// Hm_z_tab,dHm_z_tab : table of values of H-(l,eta,z) and H-'(l,eta,z).
// il_turn : l_turn = l_deb + il_turn if there is a turning l point for fixed z (see above). If there is none, l_turn = l_deb by definition.
// Gz_tab,dGz_tab : tables of G,G' values for a fixed z and varying l.
// omega : +1 or -1 so |H[omega,l_turn]| <= |H[-omega,l_turn]|. This insures H[omega,l_turn] to be minimal if F is not.
// Hz_omega_tab,dHz_omega_tab : tables of H+,H+' values for fixed z and varying l if omega = 1
//                              tables of H-,H-' values for fixed z and varying l if omega = -1
// Hz_m_omega_tab,dHz_m_omega_tab : tables of H-,H-' values for fixed z and varying l if omega = 1
//                                  tables of H+,H+' values for fixed z and varying l if omega = -1
// Fl,dFl : values of F,F' for fixed z and l. 
// H_omega_l,dH_omega_l,H_omega_lp1: values of H[omega],H'[omega] for fixed z and l, and fixed z and l+1 respectively.
// dH_omega_lp1,H_omega_lm1,dH_omega_lm1 : values of H[-omega],H'[-omega] for fixed z and l, and fixed z and l+1 respectively.                  
// Slp1,Rlp1 : Slp1 = S[l+1] = (l+1)/z + eta/(l+1) and Rlp+1 = R[l+1] = R_tab[il+1], for l = l_deb+il and il in [0:Nl-2].
// Iomega, G_factor, H_factor : i.omega and factors so H[-omega] = H[omega] - H_factor.F and G = H[omega] - G_factor.F .
//                              G_factor is i.omega for standard normalization and i.omega.(C(l,eta)^2) for non-standard normalization.
//                              H_factor is 2.i.omega for standard normalization and 2.i.omega.(C(l,eta)^2) for non-standard normalization.
// G_factor_turn : same as G_factor for l = l_turn.

void cwf_l_tables_rec_rel_helper (const int Nl,const complex<double> R_tab[],class Coulomb_wave_functions *const cwf_l_tab[],
				  const complex<double> &z,const int il_turn,
				  const complex<double> Fz_tab[],const complex<double> dFz_tab[],
				  complex<double> Hp_z_tab[],complex<double> dHp_z_tab[],
				  complex<double> Hm_z_tab[],complex<double> dHm_z_tab[],
				  complex<double> Gz_tab[],complex<double> dGz_tab[])
{
  const bool is_it_normalized = cwf_l_tab[0]->is_it_normalized;
  const complex<double> l_deb = cwf_l_tab[0]->l, eta = cwf_l_tab[0]->eta;

  cwf_l_tab[il_turn]->F_dF_init (z,Fz_tab[il_turn],dFz_tab[il_turn]);
  cwf_l_tab[il_turn]->H_dH (1,z,Hp_z_tab[il_turn],dHp_z_tab[il_turn]);
  cwf_l_tab[il_turn]->H_dH (-1,z,Hm_z_tab[il_turn],dHm_z_tab[il_turn]);

  const int omega = (abs (Hp_z_tab[il_turn]) <= abs (Hm_z_tab[il_turn])) ? (1) : (-1);
  complex<double> *const Hz_omega_tab = (omega == 1) ? (Hp_z_tab) : (Hm_z_tab),*const dHz_omega_tab = (omega == 1) ? (dHp_z_tab) : (dHm_z_tab);
  complex<double> *const Hz_m_omega_tab = (omega == 1) ? (Hm_z_tab) : (Hp_z_tab),*const dHz_m_omega_tab = (omega == 1) ? (dHm_z_tab) : (dHp_z_tab);
  
  const complex<double> I_omega(0,omega),G_factor_turn = (!is_it_normalized) ? (I_omega*exp (2.0*log_Cl_eta_calc (l_deb + il_turn,eta))) : (I_omega);
  Gz_tab[il_turn] = Hz_omega_tab[il_turn] - G_factor_turn*Fz_tab[il_turn], dGz_tab[il_turn] = dHz_omega_tab[il_turn] - G_factor_turn*dFz_tab[il_turn];

  for (int il = il_turn-1 ; il >=0 ; il--)
  {
    const complex<double> l = l_deb+il,Slp1 = (l+1)/z + eta/(l+1),Rlp1 = R_tab[il+1],H_omega_lp1 = Hz_omega_tab[il+1],dH_omega_lp1 = dHz_omega_tab[il+1];
    complex<double> &H_omega_l = Hz_omega_tab[il],&dH_omega_l = dHz_omega_tab[il];
    H_omega_l = (is_it_normalized) ? ((dH_omega_lp1 + Slp1*H_omega_lp1)/Rlp1) : ((dH_omega_lp1 + Slp1*H_omega_lp1)*(2*l+3)/(Rlp1*Rlp1));
    dH_omega_l = (is_it_normalized) ? (Slp1*H_omega_l - Rlp1*H_omega_lp1) : (Slp1*H_omega_l - (2*l + 3)*H_omega_lp1);

    const complex<double> G_factor = (!is_it_normalized) ? (I_omega*exp (2.0*log_Cl_eta_calc (l,eta))) : (I_omega),H_factor = 2*G_factor;
    const complex<double> Fl = Fz_tab[il],dFl = dFz_tab[il];
    Hz_m_omega_tab[il] = H_omega_l - Fl*H_factor, dHz_m_omega_tab[il] = dH_omega_l - dFl*H_factor;
    Gz_tab[il] = H_omega_l - G_factor*Fl, dGz_tab[il] = dH_omega_l - G_factor*dFl;
  }

  for (int il = il_turn+1 ; il < Nl ; il++)
  {
    const complex<double> l = l_deb+il,Sl = l/z + eta/l,Rl = R_tab[il],H_omega_lm1 = Hz_omega_tab[il-1],dH_omega_lm1 = dHz_omega_tab[il-1];
    complex<double> &H_omega_l = Hz_omega_tab[il],&dH_omega_l = dHz_omega_tab[il];
    H_omega_l = (is_it_normalized) ? ((Sl*H_omega_lm1 - dH_omega_lm1)/Rl) : ((Sl*H_omega_lm1 - dH_omega_lm1)/(2*l+1));
    dH_omega_l = (is_it_normalized) ? (Rl*H_omega_lm1 - Sl*H_omega_l) : (((Rl*Rl)/(2*l+1))*H_omega_lm1 - Sl*H_omega_l);

    const complex<double> G_factor = (!is_it_normalized) ? (I_omega*exp (2.0*log_Cl_eta_calc (l,eta))) : (I_omega),H_factor = 2*G_factor;
    const complex<double> Fl = Fz_tab[il],dFl = dFz_tab[il];
    Hz_m_omega_tab[il] = H_omega_l - Fl*H_factor, dHz_m_omega_tab[il] = dH_omega_l - dFl*H_factor;
    Gz_tab[il] = H_omega_l - G_factor*Fl, dGz_tab[il] = dH_omega_l - G_factor*dFl;
  }
}



// Routine calculating all F,G,H+,H- and derivative tables with recurrence relations for a table of z variables.
// -------------------------------------------------------------------------------------------------------------
// 
// It is first checked that 1+l_deb+/-i.eta is not a negative integer, as this case is not possible here.
// F,G,H+,H- and derivative are calculated with previous routines and they are checked with the wronskian between F and H[omega]
// for each l and z, with omega=+1 or -1 so |H[omega,l]| <= |H[-omega,l]|.
// 
// Variables
// ---------
// l_deb,l_end : first and last considered angular momenta.
// Nl,l : number of angular momenta, angular momentum l = l_deb+i, i in [0:Nl-1]. l_end is l_deb + Nl - 1.
// eta : Sommerfeld parameter
// Ieta,l_deb_p1_p_Ieta,l_deb_p1_m_Ieta : i.eta, 1+l+i.eta, 1+l-i.eta
// is_it_normalized : true if one uses standard normalization, false if not.
// Nz, z_tab : number of z variables, table of z variables.
// R_tab : table of R[l] values : R[l] = (2l+1).C(l,eta)/C(l-1,eta), l in [l_deb+1,l_end].
// cwf_l_tab : table of pointers to Coulomb_wave_functions classes of respective parameters l=l_deb+il, eta, is_it_normalized for il in [0:Nl-1].
// F_tab,dF_tab : table of values of F(l,eta,z) and F'(l,eta,z). F(l,eta,z) = F_tab[iz][il], F'(l,eta,z) = dF_tab[iz][il], with z = z_tab[iz] and l = l_deb+il.
// G_tab,dG_tab : table of values of G(l,eta,z) and G'(l,eta,z). G(l,eta,z) = G_tab[iz][il], G'(l,eta,z) = dG_tab[iz][il].
// Hp_tab,dHp_tab : table of values of H+(l,eta,z) and H+'(l,eta,z). H+(l,eta,z) = Hp_tab[iz][il], H+'(l,eta,z) = dHp_tab[iz][il].
// Hm_tab,dHm_tab : table of values of H-(l,eta,z) and H-'(l,eta,z). H-(l,eta,z) = Hm_tab[iz][il], H-'(l,eta,z) = dHm_tab[iz][il].
// il_turn_tab : l_turn = l_deb + il_turn_tab[iz] if there is a turning l point for a given z (see above). If there is none, l_turn = l_deb by definition.
// il_turn : l_turn = l_deb + il_turn if there is a turning l point for fixed z (see above). If there is none, l_turn = l_deb by definition.
// Fz_tab,dFz_tab : tables of F,F' values for a fixed z and varying l.

void cwf_l_tables_recurrence_relations (const complex<double> &l_deb,const int Nl,const complex<double> &eta,const bool is_it_normalized,
					const int Nz,const complex<double> z_tab[],
					complex<double> *const *const F_tab,complex<double> *const *const dF_tab,
					complex<double> *const *const G_tab,complex<double> *const *const dG_tab,
					complex<double> *const *const Hp_tab,complex<double> *const *const dHp_tab,
					complex<double> *const *const Hm_tab,complex<double> *const *const dHm_tab)
{  
  const complex<double> Ieta(-imag (eta),real (eta)), l_deb_p1_p_Ieta = 1+l_deb+Ieta, l_deb_p1_m_Ieta = 1+l_deb-Ieta;
  if ((l_deb == rint (real (l_deb))) && (real (l_deb) < 0)) cout<<"l_deb cannot be a negative integer."<<endl, exit (1);
  if ((l_deb_p1_p_Ieta == rint (real (l_deb_p1_p_Ieta))) && (real (l_deb_p1_p_Ieta) < 0)) cout<<"1+l_deb+i.eta cannot be a negative integer."<<endl, exit (1);
  if ((l_deb_p1_m_Ieta == rint (real (l_deb_p1_m_Ieta))) && (real (l_deb_p1_m_Ieta) < 0)) cout<<"1+l_deb-i.eta cannot be a negative integer."<<endl, exit (1);

  complex<double> *const R_tab = new complex<double> [Nl];
  class Coulomb_wave_functions **cwf_l_tab = new class Coulomb_wave_functions * [Nl];

  for (int il = 0 ; il < Nl ; il++)
  {
    const complex<double> l = l_deb+il;
    if (il > 0) R_tab[il] = (2*l + 1)*exp (log_Cl_eta_calc (l,eta) - log_Cl_eta_calc (l-1,eta));
    cwf_l_tab[il] = new class Coulomb_wave_functions (is_it_normalized,l,eta);
  }

  int *const il_turn_tab = new int [Nz];
  F_dF_l_tables_rec_rel_helper (Nl,Nz,z_tab,R_tab,*(cwf_l_tab[0]),*(cwf_l_tab[Nl-1]),F_tab,dF_tab,il_turn_tab);

  for (int iz = 0 ; iz < Nz ; iz++)
  { 
    const int il_turn = (il_turn_tab[iz] == -1) ? (0) : (il_turn_tab[iz]);
    const complex<double> *const Fz_tab = F_tab[iz],*const dFz_tab = dF_tab[iz];
    cwf_l_tables_rec_rel_helper (Nl,R_tab,cwf_l_tab,z_tab[iz],il_turn,Fz_tab,dFz_tab,Hp_tab[iz],dHp_tab[iz],Hm_tab[iz],dHm_tab[iz],G_tab[iz],dG_tab[iz]);
  }

  for (int il = 0 ; il < Nl ; il++) delete cwf_l_tab[il];
  delete [] cwf_l_tab;
  delete [] il_turn_tab;
  delete [] R_tab;
}




// Routine calculating all F,G,H+,H- and derivative tables with recurrence relations for a single z.
// -------------------------------------------------------------------------------------------------
// 
// It is first checked that 1+l_deb+/-i.eta is not a negative integer, as this case is not possible here.
// F,G,H+,H- and derivative are calculated with previous routines and they are checked with the wronskian between F and H[omega]
// for each l, with omega=+1 or -1 so |H[omega,l]| <= |H[-omega,l]|. 
//
// Variables
// ---------
// l_deb,l_end : first and last considered angular momenta.
// Nl,l : number of angular momenta, angular momentum l = l_deb+i, i in [0:Nl-1]. l_end is l_deb + Nl - 1.
// eta : Sommerfeld parameter
// Ieta,l_deb_p1_p_Ieta,l_deb_p1_m_Ieta : i.eta, 1+l+i.eta, 1+l-i.eta
// is_it_normalized : true if one uses standard normalization, false if not.
// z : Coulomb wave function variable
// R_tab : table of R[l] values : R[l] = (2l+1).C(l,eta)/C(l-1,eta), l in [l_deb+1,l_end].
// cwf_l_tab : table of pointers to Coulomb_wave_functions classes of respective parameters l=l_deb+il, eta, is_it_normalized for il in [0:Nl-1].
// F_tab,dF_tab : table of values of F(l,eta,z) and F'(l,eta,z). F(l,eta,z) = F_tab[il], F'(l,eta,z) = dF_tab[il], with l = l_deb+il.
// G_tab,dG_tab : table of values of G(l,eta,z) and G'(l,eta,z). G(l,eta,z) = G_tab[il], G'(l,eta,z) = dG_tab[il].
// Hp_tab,dHp_tab : table of values of H+(l,eta,z) and H+'(l,eta,z). H+(l,eta,z) = Hp_tab[il], H+'(l,eta,z) = dHp_tab[il].
// Hm_tab,dHm_tab : table of values of H-(l,eta,z) and H-'(l,eta,z). H-(l,eta,z) = Hm_tab[il], H-'(l,eta,z) = dHm_tab[il].

void cwf_l_tables_recurrence_relations (const complex<double> &l_deb,const int Nl,const complex<double> &eta,const bool is_it_normalized,
					const complex<double> z,
					complex<double> *const F_tab,complex<double> *const dF_tab,
					complex<double> *const G_tab,complex<double> *const dG_tab,
					complex<double> *const Hp_tab,complex<double> *const dHp_tab,
					complex<double> *const Hm_tab,complex<double> *const dHm_tab)
{
  const complex<double> Ieta(-imag (eta),real (eta)), l_deb_p1_p_Ieta = 1+l_deb+Ieta, l_deb_p1_m_Ieta = 1+l_deb-Ieta;
  if ((l_deb == rint (real (l_deb))) && (real (l_deb) < 0)) cout<<"l_deb cannot be a negative integer."<<endl, exit (1);
  if ((l_deb_p1_p_Ieta == rint (real (l_deb_p1_p_Ieta))) && (real (l_deb_p1_p_Ieta) < 0)) cout<<"1+l_deb+i.eta cannot be a negative integer."<<endl, exit (1);
  if ((l_deb_p1_m_Ieta == rint (real (l_deb_p1_m_Ieta))) && (real (l_deb_p1_m_Ieta) < 0)) cout<<"1+l_deb-i.eta cannot be a negative integer."<<endl, exit (1);

  cwf_l_tables_recurrence_relations (l_deb,Nl,eta,is_it_normalized,1,&z,&F_tab,&dF_tab,&G_tab,&dG_tab,&Hp_tab,&dHp_tab,&Hm_tab,&dHm_tab);
}
