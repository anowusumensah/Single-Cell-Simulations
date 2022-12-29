#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:13:09 2022

@author: cce2022
"""
import numpy as np
#from scipy.integrate import solve_ivp


def Aslanidi(t,y,stimtimes):
#def Aslanidi(t,y):    
    #global stimtimes
    #===============================================================================================
    #State variables
    V =  y[0]
    m =  y[1]
    h =  y[2]
    j =  y[3]
    mL = y[4]
    hL = y[5]
    xr = y[6]
    x  = y[7]
    y1 = y[8]
    y2 = y[9]
    xs = y[10]
    d  = y[11]
    f  = y[12]
    b  = y[13]
    g  = y[14]
    Cai= y[15]   #Ca_buffer_Cai (umolar)
    Y0 = y[16]   #Ca_buffer_Ca_Calsequestrin (millimolar)
    Y1 = y[17]   #Ca_buffer_Ca_SL (millimolar)
    Y2 = y[18]   #Ca_buffer_Ca_SLB_SL (millimolar)
    Y3 = y[19]   #Ca_buffer_Ca_SLB_jct (millimolar)
    Y4 = y[20]   #Ca_buffer_Ca_SLHigh_SL (millimolar)
    Y5 = y[21]   #Ca_buffer_Ca_SLHigh_jct (millimolar)
    Y6 = y[22]   #Ca_buffer_Ca_SR (millimolar)
    Y7 = y[23]   #Ca_buffer_Ca_jct (millimolar )
    Y11 = y[24]
    Y23 = y[25]  #Jrel_SR_I (dimensionless)
    Y24 = y[26]  #Jrel_SR_O (dimensionless)
    Y25 = y[27]  #Jrel_SR_R (dimensionless)
    Y32 = y[28]  #cytosolic_Ca_buffer_Ca_Calmodulin (millimolar)
    Y33 = y[29]  #cytosolic_Ca_buffer_Ca_Myosin (millimolar)
    Y34 = y[30]  #cytosolic_Ca_buffer_Ca_SRB (millimolar)
    Y35 = y[31]  #cytosolic_Ca_buffer_Ca_TroponinC (millimolar)
    Y36 = y[32]  #cytosolic_Ca_buffer_Ca_TroponinC_Ca_Mg (millimolar)
    Y37 = y[33]  #cytosolic_Ca_buffer_MgMyosin (millimolar)
    Y38 = y[34]  #cytosolic_Ca_buffer_MgTroponinC_Ca_Mg (millimolar)
    i =   y[35]
    #=========================================================================================================================
    # Parameters
    GNa   = 20;                      #units(mS / cm^2);
    GNaL  = 0.0162;                  #units(mS / cm^2);
    GKr   = 0.015583333333333;       #units(mS / cm^2);
    GKs   = 0.011738183745949;       #units(mS / cm^2);
    GK1   = 0.5;                     #units(mS / cm^2);
    Gto   = 0.112;                   #units(mS / cm^2);
    GCaL  = 0.27;                    #units(mS / cm^2);
    GCaT  = 0.2;                     #units(mS / cm^2);
    GCl   = 0.3;                     #units(mS / cm^2);
    GClB  = 0.0009;                  #units(mS / cm^2);
    GCaB  = 3.51820e-04;             #units(mS / cm^2);
    GNaB  = 2.97000e-05;             #units(mS / cm^2);
    GKB   = 5.00000e-05;             #units(mS / cm^2);
    Nae     = 140.0;                 #units(mM);
    Nai     = 8.8;                   #units(mM);
    Ki      = 135.0;                 #units(mM);
    Ke      = 5.4;                   #units(mM);
    Cle     = 150.0;                 #units(mM);
    Cli     = 30.0;                  #units(mM);
    Mgi     = 1.0;                   #units(mM);
    Cae     = 1.8;                   #units(mM);
    T = 310.0;                       #units(K);
    INaCamax = 4.5;
    kmcaact  = 0.000125;  #units(mM);   
    kmnai1   = 12.3;      #units(mM);
    kmnao    = 87.5;      #units(mM);
    kmcai    = 0.0036;    #units(mM);       
    kmcao    = 1.3;       #units(mM);
    nu       = 0.35;      #position of energy barrier for INaCa
    ksat     = 0.27;      #saturation factor for INaCa at negatiVole potentials
    PNaK    = 0.01833;
    Ca_handling = 1;      #by default, Ca2+ handling is turned on
    ICaPHill = 1.6;
    ICaPVmf  = 0.0538 * 1.25;
    ICaPKmf  = 0.5;
    #============================================================================================================
    F     = 96486.7;        #units(J / kmol / K);
    R     = 8314.3;         #units(C / mol);
    RTonF = R * T / F;          #units(mV);
    PI    = 3.14159265359;  #unitless  
    Ca_buffer_Bmax_Calsequestrin  = 0.14;       #mM
    Ca_buffer_Bmax_SLB_SL         = 0.0374;     #mM
    Ca_buffer_Bmax_SLB_jct        = 0.0046;     #mM
    Ca_buffer_Bmax_SLHigh_SL      = 0.0134;     #mM
    Ca_buffer_Bmax_SLHigh_jct     = 0.00165;    #mM
    Ca_buffer_koff_Calsequestrin  = 65.0;       #ms^-1
    Ca_buffer_koff_SLB            = 1.3;        #ms^-1
    Ca_buffer_koff_SLHigh         = 30.0e-3;    #ms^-1
    Ca_buffer_kon_Calsequestrin   = 100.0;      #mM^-1 ms^-1
    Ca_buffer_kon_SL              = 100.0;      #mM^-1 ms^-1
    JleaKSR_KSRleak               = 5.348e-6;   #per_millisecond
    JpumPSR_H                     = 1.787;      #dimensionless
    JpumPSR_Kmf                   = 0.000246;   #millimolar
    JpumPSR_Kmr                   = 1.7;        #millimolar
    #JpumPSR_Q10_SRCaP             = 2.6;        #dimensionless
    JpumPSR_V_max                 = 286.0e-6;   #millimolar_per_millisecond
    Jrel_SR_EC50_SR               = 0.45;       #millimolar
    Jrel_SR_HSR                   = 2.5;        #dimensionless
    Jrel_SR_Max_SR                = 15.0;       #dimensionless
    Jrel_SR_Min_SR                = 1.0;        #dimensionless
    Jrel_SR_kiCa                  = 0.5;        #per_millimolar_per_millisecond
    Jrel_SR_kim                   = 0.005;      #per_millisecond
    Jrel_SR_koCa                  = 10.0;       #per_millimolar2_per_millisecond
    Jrel_SR_kom                   = 0.06;       #per_millisecond
    Jrel_SR_ks                    = 25.0;       #per_millisecond
    cytosolic_Ca_buffer_Bmax_Calmodulin         = 0.024;      #millimolar
    cytosolic_Ca_buffer_Bmax_Myosin_Ca          = 0.14;       #millimolar
    cytosolic_Ca_buffer_Bmax_Myosin_Mg          = 0.14;       #millimolar
    cytosolic_Ca_buffer_Bmax_SRB                = 0.0171;     #millimolar
    cytosolic_Ca_buffer_Bmax_TroponinC          = 0.07;       #millimolar
    cytosolic_Ca_buffer_Bmax_TroponinC_Ca_MGCa  = 0.14;       #millimolar
    cytosolic_Ca_buffer_Bmax_TroponinC_Ca_MgMg  = 0.14;       #millimolar
    cytosolic_Ca_buffer_koff_Calmodulin         = 238.0e-3;   #per_millisecond
    cytosolic_Ca_buffer_koff_Myosin_Ca          = 0.46e-3;    #per_millisecond
    cytosolic_Ca_buffer_koff_Myosin_Mg          = 0.057e-3;   #per_millisecond
    cytosolic_Ca_buffer_koff_SRB                = 60.0e-3;    #per_millisecond
    cytosolic_Ca_buffer_koff_TroponinC          = 19.6e-3;    #per_millisecond
    cytosolic_Ca_buffer_koff_TroponinC_Ca_MGCa  = 0.032e-3;   #per_millisecond
    cytosolic_Ca_buffer_koff_TroponinC_Ca_MgMg  = 3.33e-3;    #per_millisecond
    cytosolic_Ca_buffer_kon_Calmodulin          = 34.0;       #per_millimolar_per_millisecond
    cytosolic_Ca_buffer_kon_Myosin_Ca           = 13.8;       #per_millimolar_per_millisecond
    cytosolic_Ca_buffer_kon_Myosin_Mg           = 15.7e-3;    #per_millimolar_per_millisecond
    cytosolic_Ca_buffer_kon_SRB                 = 100.0;      #per_millimolar_per_millisecond
    cytosolic_Ca_buffer_kon_TroponinC           = 32.7;       #per_millimolar_per_millisecond
    cytosolic_Ca_buffer_kon_TroponinC_Ca_MGCa   = 2.37;       #per_millimolar_per_millisecond
    cytosolic_Ca_buffer_kon_TroponinC_Ca_MgMg   = 3.0e-3;     #per_millimolar_per_millisecond 
# =============================================================================
#     ion_diffusion_A_SL_cytosol     = 1.3e-4;     #cm2
#     ion_diffusion_A_jct_SL         = 3.01e-6;    #cm2
#     ion_diffusion_D_Ca_SL_cytosol  = 1.22e-6;    #dm2_per_second
#     ion_diffusion_D_Ca_jct_SL      = 1.64e-6;    #dm2_per_second
#     ion_diffusion_D_NaSL_cytosol   = 1.79e-5;    #dm2_per_second
#     ion_diffusion_D_Najct_SL       = 1.09e-5;    #dm2_per_second
#     ion_diffusion_x_SL_cytosol     = 0.45;       #micrometre
#     ion_diffusion_x_jct_SL         = 0.5;        #micrometre
# =============================================================================
    Cm_per_area   = 2.0e-6;     #farad_per_cm2
    cell_length   = 100.0;      #micrometre
    cell_radius   = 10.25;      #micrometre
    Vol_Cell      = PI * pow(cell_radius / 1000.0,2) * cell_length / pow(1000.0,3);         #units(L);
    Vol_cytosol   = 0.65 * Vol_Cell;                                                        #units(L);
    Vol_SR        = 0.035 * Vol_Cell;                                                       #units(L);
    Vol_SL        = 0.02 * Vol_Cell;                                                        #units(L);
    Vol_jct       = 0.00051 * Vol_Cell;                                                     #units(L);
    Cm            = Cm_per_area * 2.0 * cell_radius / 10000.0 * PI * cell_length / 10000.0; #units(F);
    #=================================================================================================================================================================================================
    # Simple stimulus current
    Istim = 0;
    stimduration=1.0;
    for count in np.arange(np.shape(stimtimes)[0]):
        if (stimtimes[count] <= t < stimtimes[count]+stimduration):
            Istim = -60;
    Cai_mM = Cai / 1000.;
    
    #Nernst Equilibrium Potential
    E_Na = RTonF * np.log(Nae / Nai);                                #units(mV);
    E_K  = RTonF * np.log(Ke  / Ki );                                #units(mV);
    E_Ks = RTonF * np.log((Ke + PNaK * Nae)/(Ki + PNaK * Nai));      #units(mV);
    
    #=====================================================================================================================================
    #Ionic currents
    #=====================================================================================================================================
    #Fast INa
    a_m = 0.32 * (V+47.13)/(1- np.exp (-0.1* (V+47.13)));
    b_m = 0.08 * np.exp (-V / 11);
    if (V >= -40): 
        a_h =  0; 
        b_h = 1/(0.13 * (1+np.exp ((V+10.66)/-11.1)));
        a_j = 0;
        b_j =  0.3 * np.exp (-2.535 * 10e-7 * V)/(1+ np.exp (-0.1* (V+32)));
    else:
        a_h = 0.135 * np.exp ((80+V)/-6.8);
        b_h = 0.3 * np.exp (-2.535 * 10e-7 * V)/(1+ np.exp (-0.1* (V+32)));
        a_j = ((-1.2714 * 10e5 * np.exp (0.2444* V) - (3.474 * 10e-5 * np.exp (-0.0439 * V)) * (V+37.78)))/(1+np.exp (0.311 * (V+79.23)));
        b_j = 0.1212 * np.exp (-0.01052* V)/(1+np.exp (-0.1378 * (V+40.14)));
        
    m_inf = a_m/(a_m + b_m);
    tau_m = 1./(a_m + b_m);
    h_inf = a_h/(a_h + b_h);
    tau_h = 1./(a_h + b_h);
    j_inf = a_j/(a_j + b_j);
    tau_j = 1./(a_j + b_j);
    diff_m = (m_inf - m) /tau_m;
    diff_h = (h_inf - h) /tau_h;
    diff_j = (j_inf - j) /tau_j;
    INa = GNa * m * m * m * h * j * (V - E_Na);
    
    #Slow INa
    if (abs(V+47.13) < 0.05):
        a_mL = 0.32 / 0.1;
    else:
        a_mL = 0.32*(V+47.13)/(1-np.exp(-0.1*(V+47.13)));
    b_mL   = 0.08 * np.exp(-V / 11);  
    mL_inf = a_mL/(a_mL + b_mL);
    tau_mL = 1./(a_mL + b_mL);
    tau_hL = 132.4+ 112.8 * np.exp(0.02325 * V);
    hL_inf = 1/(1+np.exp((V-22.0+91.0)/6.1));
    diff_mL = (mL_inf - mL) /tau_mL;
    diff_hL = (hL_inf - hL) /tau_hL;
    INaL  = GNaL * mL * hL * (V - E_Na);

    # Background INa current
    INaB  = GNaB * (V - E_Na);
    
    #INaK
    sigma = (np.exp(Nae / 67.3) - 1.0)/7.0;
    fNaK  = 1.0/(1.0 + 0.1245 * np.exp(-0.1 * V/(RTonF)) + 0.0365 * sigma * np.exp(-1.0 * V/(RTonF)));
    INaK = 0.61875 * fNaK * 1/(1+pow(10 / Nai,2))*(Ke/(Ke+1.5));
    
    #INCX
    Nai3 = Nai * Nai * Nai;
    Nae3 = Nae * Nae * Nae;
    kmnai13     = kmnai1 * kmnai1 * kmnai1;
    kmnao3      = kmnao * kmnao * kmnao;
    exPnuVFRT   = np.exp(nu * V * F/(R * T));
    allo        = 1/(1+pow((kmcaact/(1.5 * Cai_mM)),2));
    exPnum1VFRT = np.exp((nu-1)*V * F/(R * T));
    num         = 0.4 * INaCamax*(Nai3 * Cae * exPnuVFRT-Nae3 * 1.5 * Cai_mM * exPnum1VFRT);
    denommult   = 1+ksat * np.exp((nu-1)*V * F/(R * T));
    denomterm1  = kmcao * Nai3+kmnao3 * 1.5 * Cai_mM+kmnai13 * Cae*(1+1.5 * Cai_mM / kmcai);
    denomterm2  = kmcai * Nae3*(1+Nai3 / kmnai13)+Nai3 * Cae+Nae3 * 1.5 * Cai_mM;
    deltaE      = num/(denommult*(denomterm1+denomterm2));
    INCX = allo * deltaE;

    #IKr
    r_infty = 1.0/(1.0 + np.exp(V / 50));
    xr_inf = 1.0/(1.0 + np.exp(-1.0 * (V + 20.0)/10.5));
# =============================================================================
#     Actual formulation from the paper
#     f1 = 1.0*(1- np.exp(-0.1233*(V + 7)))/0.00138*(V + 7.0)
#     f2 = 0.00061*(V + 10.0)/ (np.exp(0.145*(V + 10.0)) - 1.0)
#     tau_xr = f1 + f2
#     
# =============================================================================
    # tau_xr formulation in openCARP
    if (abs(V+7)  < 0.05):
        tau_xr1 = 0.00138 / 0.123
    else:
        tau_xr1 = 0.00138 * (V + 7.0)/(1.0 - np.exp(-0.123 * (V + 7.0)))
        
    if (abs(V+10) < 0.05):
        tau_xr2 = 0.00061 / 0.145
    else:
        tau_xr2 = 0.00061 * (V + 10.0)/(np.exp(0.145 * (V + 10.0)) - 1.0)    
    tau_xr = 1.0/(tau_xr1 + tau_xr2)
    diff_xr = (xr_inf - xr)/tau_xr;
    IKr = GKr * xr * r_infty * (V - E_K);
    
    #IKs
    xs_inf = 1.1/(1.0+np.exp(-1.0*(V-1.5)/20));
    tau_xs = (600.0/(1.0+np.exp((V-20)/15.0)) + 250.0);
    diff_xs = (xs_inf - xs)/tau_xs;
    IKs = GKs * xs * (V - E_Ks);
    
    #IK1
    aa_k1 = 0.3/(1.0 + np.exp(0.2385 * (V - E_K - 59.215)));
    bb_k1 = (0.49124 * np.exp(0.08032 * (V - E_K + 5.476)) + np.exp(0.06175 * (V - E_K - 594.31)))/(1.0 + np.exp(-0.5143 * (V - E_K + 4.753)));
    k1_infty = aa_k1/(aa_k1 + bb_k1);
    IK1 = GK1 * (k1_infty + 0.008) * (V - E_K);
    
    #Plateau K+ current
    IKPKp = 1.0/(1.0 + np.exp((7.488 - V)/5.98)); 
    IKp    = 0.001 * IKPKp * (V - E_K);
    
    #Background Potassium current
    IKB  = GKB * (V - E_K);
    
    #Ito
    aa_x  = (0.0451 * np.exp (0.85 * 0.03577 * V));
    bb_x  = (0.0989 * np.exp (-0.85 * 0.0623 * V));
    x_inf = aa_x / (aa_x + bb_x);
    tau_x = 0.2  / (aa_x + bb_x);
    aa_y   = (0.05415 * np.exp (-(V+12.5)/15)/(1+0.051335 * np.exp(-(V+12.5)/15)));
    bb_y   = (0.05415 * np.exp ((V+33.5)/15)/(1+0.051335 * np.exp((V+33.5)/15)));
    y1_inf = aa_y/(aa_y+bb_y);
    y2_inf = y1_inf;
    tau_y1 = 0.7*(15+20.0/(aa_y+bb_y));
    tau_y2 = 4.0/(aa_y+bb_y);
    diff_x = (x_inf - x)/tau_x;
    diff_y1 = (y1_inf-y1)/tau_y1;
    diff_y2 = (y2_inf -y2)/tau_y2;
    Ito = Gto * x * (y1 * 0.75 + y2 * 0.25) * (V - E_K);
    
    #ICaL
    E_CaL = 60.0;   #units(mV);
    d_inf = (1/(1.0+np.exp(-(V-4.0)/6.74)));
    tau_d = (0.59+0.8 * np.exp(0.052*(V+13))/(1+np.exp(0.132*(V+13))));
    f_inf = 1./(1.0 + np.exp((V+25.)/10.));
    tau_f = 0.005 * np.square((V-2.5))+4.0;
    diff_d = (d_inf-d)/tau_d;
    diff_f = (f_inf-f)/tau_f;
    ICaL = GCaL * d * f * (1-Y11) * (V - E_CaL);

    #ICaT
    E_CaT = 50.0;   #units(mV);
    b_inf = 1/(1+ np.exp (-(V+28)/6.1));
    aa_b = 1.068 * np.exp((V+16.3)/30.0);
    bb_b  = 1.068 * np.exp(-(V+16.3)/30.0);
    tau_b = 1.0/(aa_b+bb_b);
    g_inf = 1/(1+np.exp((V+60)/6.6));
    aa_g  = 0.015 * np.exp(-(V+71.7)/83.3);
    bb_g  = 0.015 * np.exp((V+71.7)/15.4); 
    tau_g = 1.0/(aa_g+bb_g);
    diff_b = (b_inf-b)/tau_b;
    diff_g = (g_inf-g)/tau_g;
    ICaT = GCaT * b * g * (V - E_CaT);
    
    #ICa = ICaT + ICaL
    ICa = ICaL + ICaT;
    
    #Background ICa
    E_Ca  = (RTonF)/2.0 * np.log(Cae/(Cai_mM));
    ICaB = GCaB * (V - E_Ca);
    
    #Ca pump current
    ICaP = 0.5 * ICaPVmf / (1.0 + pow((ICaPKmf/(Cai_mM * 1000)), ICaPHill));
    
    #Ca activated Cl- current
    E_Cl = RTonF * np.log(Cli / Cle);
    a_inf = 1.0/(1.0 + np.exp(-(V + 5.0)/10.0));
    i_inf = 1.0/(1.0 + np.exp((V + 75.0)/10.0));
    tau_i = (10.0 / (1.0 + np.exp((V + 33.5)/10.0)) + 10.0)*2;
    diff_i = (i_inf-i)/tau_i;
    ICl  = GCl * a_inf * i / (1.0 + 0.1/(Cai_mM * 1000)) * (V - E_Cl);
    
    #Background chloride current
    IClB = GClB* (V - E_Cl);
    
    #Iion
    Iion = INa + INaL + IKr + IKs + IK1 + IKp + Ito + ICaL + ICaT + ICl+ IClB + ICaB + INaB + IKB + INaK + ICaP + INCX;
    Itotal = Iion + Istim;
    diff_V =-Itotal;

    
    #Calcium Handling
    Ca_buffer_dCalsequestrin = Ca_buffer_kon_Calsequestrin * Y6 * (Ca_buffer_Bmax_Calsequestrin * Vol_cytosol / Vol_SR - Y0) - Ca_buffer_koff_Calsequestrin * Y0;
    diff_Y0  = Ca_buffer_dCalsequestrin;
    Ca_buffer_dCa_SLB_SL     = Ca_buffer_kon_SL * Y1 * (Ca_buffer_Bmax_SLB_SL * Vol_cytosol / Vol_SL-Y2)-Ca_buffer_koff_SLB * Y2;                      
    Ca_buffer_dCa_SLB_jct    = Ca_buffer_kon_SL * Y7 * (Ca_buffer_Bmax_SLB_jct * 0.1 * Vol_cytosol / Vol_jct-Y3)-Ca_buffer_koff_SLB * Y3;                          
    Ca_buffer_dCa_SLHigh_SL  = Ca_buffer_kon_SL * Y1 * (Ca_buffer_Bmax_SLHigh_SL * Vol_cytosol / Vol_SL-Y4) -Ca_buffer_koff_SLHigh * Y4;
    Ca_buffer_dCa_SLHigh_jct = Ca_buffer_kon_SL * Y7 * (Ca_buffer_Bmax_SLHigh_jct * 0.1 * Vol_cytosol / Vol_jct-Y5)-Ca_buffer_koff_SLHigh * Y5;
    diff_Y2 = Ca_buffer_dCa_SLB_SL;
    diff_Y3 = Ca_buffer_dCa_SLB_jct;
    diff_Y4 = Ca_buffer_dCa_SLHigh_SL;
    diff_Y5 = Ca_buffer_dCa_SLHigh_jct;
    Ca_buffer_dCa_jct_tot_bound = Ca_buffer_dCa_SLB_jct + Ca_buffer_dCa_SLHigh_jct;
    Ca_buffer_dCa_SL_tot_bound  = Ca_buffer_dCa_SLB_SL  + Ca_buffer_dCa_SLHigh_SL;
    Ca_buffer_ICa_jct_tot = ICa;
    Ca_buffer_ICaSL_tot   = -2.0 * INCX + 1.0 * ICaB + 1.0 * ICaP;
    JpumPSR_j_pumPSR = 2.0 * JpumPSR_V_max * Vol_cytosol / Vol_SR*(pow(Cai_mM / JpumPSR_Kmf, JpumPSR_H)\
                     -pow(Y6 / JpumPSR_Kmr, JpumPSR_H))/(1.0+pow(Cai_mM / JpumPSR_Kmf, JpumPSR_H) \
                     +pow(Y6 / JpumPSR_Kmr, JpumPSR_H));
    JleaKSR_j_leaKSR = 0.5 * JleaKSR_KSRleak*(Y6-Y7);
    Jrel_SR_j_rel_SR = 2.0 * Jrel_SR_ks * Y24*(Y6-Y7);
    diff_Y6  = JpumPSR_j_pumPSR-(JleaKSR_j_leaKSR * Vol_cytosol / Vol_SR+Jrel_SR_j_rel_SR)-Ca_buffer_dCalsequestrin;
    ion_diffusion_J_Ca_jct_SL = (Y7 - Y1) * 8.2413e-13;
    ion_diffusion_J_Ca_SL_cytosol = (Y1 - Cai_mM) * 3.7243e-12;
    diff_Y1 = -0.5 * Ca_buffer_ICaSL_tot * Cm / (Vol_SL * 2.0 * F)\
        + (ion_diffusion_J_Ca_jct_SL - ion_diffusion_J_Ca_SL_cytosol) / Vol_SL - 1.0 * Ca_buffer_dCa_SL_tot_bound;
    cytosolic_Ca_buffer_dCa_TroponinC = cytosolic_Ca_buffer_kon_TroponinC * Cai_mM*(cytosolic_Ca_buffer_Bmax_TroponinC-Y35)-cytosolic_Ca_buffer_koff_TroponinC * Y35;                       
    cytosolic_Ca_buffer_dCa_TroponinC_Ca_Mg = cytosolic_Ca_buffer_kon_TroponinC_Ca_MGCa * Cai_mM*(cytosolic_Ca_buffer_Bmax_TroponinC_Ca_MGCa-(Y36+Y38))-cytosolic_Ca_buffer_koff_TroponinC_Ca_MGCa * Y36;
    cytosolic_Ca_buffer_dMgTroponinC_Ca_Mg = cytosolic_Ca_buffer_kon_TroponinC_Ca_MgMg * Mgi*(cytosolic_Ca_buffer_Bmax_TroponinC_Ca_MgMg-(Y36+Y38))-cytosolic_Ca_buffer_koff_TroponinC_Ca_MgMg * Y38;
    cytosolic_Ca_buffer_dCa_Calmodulin = cytosolic_Ca_buffer_kon_Calmodulin * Cai_mM*(cytosolic_Ca_buffer_Bmax_Calmodulin-Y32)-cytosolic_Ca_buffer_koff_Calmodulin * Y32;
    cytosolic_Ca_buffer_dCa_Myosin = cytosolic_Ca_buffer_kon_Myosin_Ca * Cai_mM*(cytosolic_Ca_buffer_Bmax_Myosin_Ca-(Y33+Y37))-cytosolic_Ca_buffer_koff_Myosin_Ca * Y33;
    cytosolic_Ca_buffer_dMgMyosin = cytosolic_Ca_buffer_kon_Myosin_Mg * Mgi*(cytosolic_Ca_buffer_Bmax_Myosin_Mg-(Y33+Y37))-cytosolic_Ca_buffer_koff_Myosin_Mg * Y37;
    cytosolic_Ca_buffer_dCa_SRB = cytosolic_Ca_buffer_kon_SRB * Cai_mM*(cytosolic_Ca_buffer_Bmax_SRB-Y34)-cytosolic_Ca_buffer_koff_SRB * Y34;
    cytosolic_Ca_buffer_dCa_cytosol_tot_bound = cytosolic_Ca_buffer_dCa_TroponinC+cytosolic_Ca_buffer_dCa_TroponinC_Ca_Mg+cytosolic_Ca_buffer_dMgTroponinC_Ca_Mg+cytosolic_Ca_buffer_dCa_Calmodulin+cytosolic_Ca_buffer_dCa_Myosin+cytosolic_Ca_buffer_dMgMyosin+cytosolic_Ca_buffer_dCa_SRB;
    Jrel_SR_kCaSR = Jrel_SR_Max_SR-(Jrel_SR_Max_SR-Jrel_SR_Min_SR)/(1.0+pow(Jrel_SR_EC50_SR / Y6, Jrel_SR_HSR));
    
    diff_Y32 = cytosolic_Ca_buffer_dCa_Calmodulin;
    diff_Y33 = cytosolic_Ca_buffer_dCa_Myosin;
    diff_Y34 = cytosolic_Ca_buffer_dCa_SRB;
    diff_Y35 = cytosolic_Ca_buffer_dCa_TroponinC;
    diff_Y36 = cytosolic_Ca_buffer_dCa_TroponinC_Ca_Mg;
    diff_Y37 = cytosolic_Ca_buffer_dMgMyosin;
    diff_Y38 = cytosolic_Ca_buffer_dMgTroponinC_Ca_Mg;
        
    diff_Y7  = -0.5 * Ca_buffer_ICa_jct_tot * Cm /(Vol_jct * 2.0 * F)\
         -ion_diffusion_J_Ca_jct_SL / Vol_jct \
         +Jrel_SR_j_rel_SR * Vol_SR / Vol_jct \
         +JleaKSR_j_leaKSR * Vol_cytosol / Vol_jct \
         -1.0 * Ca_buffer_dCa_jct_tot_bound;
         
    Jrel_SR_koSRCa = Jrel_SR_koCa / Jrel_SR_kCaSR;
    Jrel_SR_kiSRCa = Jrel_SR_kiCa * Jrel_SR_kCaSR;
    Jrel_SR_RI     = 1.0-Y25-Y24-Y23;  
    
    diff_Y23 = Jrel_SR_kiSRCa * Y7 * Y24-Jrel_SR_kim * Y23-(Jrel_SR_kom * Y23-Jrel_SR_koSRCa * Y7 * Y7 * Jrel_SR_RI);
    diff_Y24 = Jrel_SR_koSRCa * Y7 * Y7 * Y25-Jrel_SR_kom * Y24-(Jrel_SR_kiSRCa * Y7 * Y24-Jrel_SR_kim * Y23);
    diff_Y25 = Jrel_SR_kim * Jrel_SR_RI-Jrel_SR_kiSRCa * Y7 * Y25-(Jrel_SR_koSRCa * Y7 * Y7 * Y25-Jrel_SR_kom * Y24);

    if ( Ca_handling == 1 ):
       diff_Cai  = \
       (-JpumPSR_j_pumPSR * Vol_SR / Vol_cytosol+ion_diffusion_J_Ca_SL_cytosol / Vol_cytosol-1.0 * cytosolic_Ca_buffer_dCa_cytosol_tot_bound)*1000.;
       diff_Y11 = 0.7 * Y7*(1.0-Y11)-11.9e-3 * Y11;
    else: 
       diff_Cai = 0;
       diff_Y11 = 0;    
    
    return diff_V,diff_m,diff_h,diff_j,diff_mL,diff_hL,diff_xr,diff_x,diff_y1,diff_y2,diff_xs,diff_d,diff_f,diff_b,diff_g,diff_Cai,diff_Y0,\
        diff_Y1,diff_Y2,diff_Y3,diff_Y4,diff_Y5,diff_Y6,diff_Y7,diff_Y11,diff_Y23,diff_Y24,diff_Y25,diff_Y32,diff_Y33,diff_Y34,diff_Y35,diff_Y36,\
            diff_Y37,diff_Y38,diff_i
            
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
