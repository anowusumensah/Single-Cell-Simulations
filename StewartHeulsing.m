function out = StewartHeulsing(t,y,stimtimes)
% Function written by Anthony Owusu-Mensah
% This function simulates a ventricular myocyte that is
% coupled to a Purkinje myocyte via gap junctional resistance 
% This simulation was to mimic the paper by Huelsing et. al (DOI:
% 10.1152/ajpheart.1998.274.4.H1163)
    
% % Initial values  
% % Ventriclular parameters
%%=============================================================================
    V =     y(1);
    y_in =  y(2);
    Xr1 =   y(3);
    Xr2 =   y(4);
    Xs =    y(5);
    m =     y(6);
    h =     y(7);
    j =     y(8);
    d =     y(9);
    f =     y(10);
    f2 =    y(11);
    fCass = y(12);
    s =     y(13);
    r =     y(14);
    Cai =   y(15);
    Ca_SR = y(16);
    Cass =  y(17);
    R_prime = y(18);
    Nai =   y(19);
    Ki =    y(20);
% %============================================================================================================
    % % PurKinje Parameters 
    Vm =     y(21);
    y_inm =  y(22);
    Xr1m =   y(23);
    Xr2m =   y(24);
    Xsm =    y(25);
    mm =     y(26);
    hm =     y(27);
    jm =     y(28);
    dm =     y(29);
    fm =     y(30);
    f2m =    y(31);
    fCassm = y(32);
    sm =     y(33);
    rm =     y(34);
    Caim =   y(35);
    Ca_SRm = y(36);
    Cassm =  y(37);
    R_primem = y(38);
    Naim =   y(39);
    Kim =    y(40);
%%===================================================================================================================================

    % Constants
    stimDur = 0.5; % Stimulus duration
% %     bcl = 1000; % basic cycle length
    ifMulv = 0;
    ifMulp = 0;
    K = 1.26e-4;
    GNa   = 130.5744;  %units(nS/pF);
    GCaL  = 3.98e-5;   %units(L/s/F);
    Gto   = 0.08184;   %units(nS/pF);
    GKr   = 0.0918;    %units(nS/pF);
    GKs   = 0.2352;    %units(nS/pF);
    Gsus  = 0.0227;    %units(nS/pF);
    Gf_Na = 0.0145654; %units(nS/pF);
    Gf_K  = 0.0234346; %units(nS/pF);
    GK1   = 0.065;     %units(nS/pF);
    maxINaCa = 1000;   %units(pA/pF);
    maxINaK  = 2.724;  %units(pA/pF);
    GpCa  = 0.1238;    %units(pA/pF);
    Gbna  = 0.00029;   %units(nS/pF);
    Gbca  = 0.000592;  %units(nS/pF);
    GpK   = 0.0146;    %units(nS/pF);

    % % Constants
    Rj = 110e6;%units(ohm); % Coupling resistance btw VM and PC
    Cm = 0.185; %units(uF);
    R = 8314.472; %units(J/K/mol);
    T = 310;    %units(K);
    F = 96485.3415; %units(C/mmol);
    Nao = 140; %units(mmol/L);
    Cao = 2.0; %units(mmol/L);
    Ko  = 5.4; %units(mmol/L);
    Buf_c = 0.2; %units(mmol/L);
    min_sr = 1; %units(unitless);
    k1_prime = 0.15; %units(L^2/mmol^2/ms);
    k4 = 0.005; %units(unitless/ms);
    V_xfer = 0.0038; %units(unitless/ms);
    K_up = 0.00025; %units(mmol/L);
    V_sr = 0.001094; %units(um^3);
    V_ss = 5.468e-5; %units(um^3);
    k3 = 0.06; %units(unitless/ms);
    Vmax_up = 0.006375; %units(mmol/L/ms);
    P_kna = 0.03; %units(unitless);
    K_pCa = 0.0005; %units(mmol/L);
    V_c = 0.016404; %units(um^3);
    V_rel = 0.102; %units(unitless/ms);
    V_leak = 0.00036; %units(unitless/ms);
    Buf_ss = 0.4; %units(mmol/L);
    Buf_sr = 10; %units(mmol/L);
    K_mNa = 40; %units(mmol/L);
    EC = 1.5; %units(mmol/L);
    alpha = 2.5; %units(unitless);
    k2_prime = 0.045; %units(L/mmol/ms);
    K_buf_sr = 0.3; %units(mmol/L);
    K_buf_ss = 0.00025; %units(mmol/L);
    K_sat = 0.1; %units(unitless);
    K_buf_c = 0.001; %units(mmol/L);
    Km_Nai = 87.5; %units(mmol/L);
    max_sr = 2.5; %units(unitless);
    Km_Ca = 1.38; %units(mmol/L);
    K_mk = 1; %units(mmol/L);
    gamma = 0.35; %units(unitless);
% % Constants
  %%=============================================================================================================================
% % Ventricular Part of the simulation
%%==============================================================================================================================  
% %     reversal_potentials
    E_K = (((R*T)/F)*log((Ko/Ki))); %units(mV);
    E_Ca = (((0.5*R*T)/F)*log((Cao/Cai))); %units(mV);
    E_Ks = (((R*T)/F)*log(((Ko+(P_kna*Nao))/(Ki+(P_kna*Nai))))); %units(mV);
    E_Na = (((R*T)/F)*log((Nao/Nai))); %units(mV);

    
% %     h_inperpolarization_activated_current
    If_Na =    ifMulv*(y_in*Gf_Na*(V - E_Na)); %units(pA/pF);
    If_K =    ifMulv*(y_in*Gf_K*(V - E_K)); %units(pA/pF);
    If =    ifMulv*(If_Na+If_K); %units(pA/pF);

% %     hyperpolarization_activated_current_y_gate
    y_inf = (1./(1.+exp(((V+80.6)/6.8)))); %units(unitless);
    aalpha_y = (1.*exp((-2.9 - (0.04*V)))); %units(unitless/ms);
    bbeta_y = (1.*exp((3.6+(0.11*V)))); %units(unitless/ms);
    tau_y = (4000./(aalpha_y+bbeta_y)); %units(ms);
    diff_y = (y_inf - y_in)/tau_y;

% % inward_rectifier_potassium_current
    xK1_inf = (1./(1.+exp((0.1*(V+75.44))))); %units(unitless);
    IK1 = (3*GK1*xK1_inf*((V - 8.) - E_K)); %units(pA/pF);
    
% % rapid_time_dependent_potassium_current
    IKr = (2*GKr*sqrt((Ko/5.4))*Xr1*Xr2*(V - E_K)); %units(pA/pF);


% %     rapid_time_dependent_potassium_current_Xr1_gate
    xr1_inf = (1./(1.+exp(((-26. - V)/7.)))); %units(unitless);
    alpha_xr1 = (450./(1.+exp(((-45. - V)/10.)))); %units(unitless);
    beta_xr1 = (6./(1.+exp(((V+30.)/11.5)))); %units(unitless);
    tau_xr1 = (1.*alpha_xr1*beta_xr1); %units(ms);
    diff_Xr1 = ((xr1_inf - Xr1)/tau_xr1); %units(unitless/ms);

% %    rapid_time_dependent_potassium_current_Xr2_gate
    alpha_xr2 = (3./(1.+exp(((-60. - V)/20.)))); %units(unitless);
    xr2_inf = (1./(1.+exp(((V+88.)/24.)))); %units(unitless);
    beta_xr2 = (1.12/(1.+exp(((V - 60.)/20.)))); %units(unitless); 
    tau_xr2 = (1.*alpha_xr2*beta_xr2); %units(ms);
    diff_Xr2 = ((xr2_inf - Xr2)/tau_xr2); %units(unitless/ms);
    
% %     slow_time_dependent_potassium_current
    IKs = (2*GKs*(Xs*Xs)*(V - E_Ks)); %units(pA/pF);

% %     slow_time_dependent_potassium_current_Xs_gate
    alpha_xs = (1400./sqrt((1.+exp(((5. - V)/6.))))); %units(unitless);
    beta_xs = (1./(1.+exp(((V - 35.)/15.)))); %units(unitless);
    tau_xs = ((1.*alpha_xs*beta_xs)+80.); %units(ms);
    xs_inf = (1./(1.+exp(((-5. - V)/14.)))); %units(unitless);
    diff_Xs = ((xs_inf - Xs)/tau_xs); %units(unitless/ms);

% %     fast_sodium_current
    INa = (0.5*GNa*(m*m*m)*h*j*(V - E_Na)); %units(pA/pF);

% % fast_sodium_current_m_gate
    bbeta_m = ((0.1/(1.+exp(((V+35.)/5.))))+(0.1/(1.+exp(((V - 50.)/200.))))); %units(unitless);
    aalpha_m = (1./(1.+exp(((-60. - V)/5.)))); %units(unitless);
    m_inf = (1./((1.+exp(((-56.86 - V)/9.03)))*(1.+exp(((-56.86 - V)/9.03))))); %units(unitless);
    tau_m = (1.*aalpha_m*bbeta_m); %units(ms);
    diff_m = (m_inf - m)/tau_m;

% %     fast_sodium_current_h_gate
     if (V<-40.) 
       bbeta_h = ((2.7*exp((0.079*V)))+(310000.*exp((0.3485*V)))); 
     else
          bbeta_h =(0.77/(0.13*(1.+exp(((V+10.66)/-11.1))))); %units(unitless/ms);
     end
     if (V<-40.)
        aalpha_h = (0.057*exp((-(V+80.)/6.8)));
     else
        aalpha_h =  0.; %units(unitless/ms;
     end
    tau_h = (1./(aalpha_h+bbeta_h)); %units(ms);
    h_inf = (1./((1.+exp(((V+71.55)/7.43)))*(1.+exp(((V+71.55)/7.43))))); %units(unitless);
    diff_h = (h_inf - h)/tau_h;
    
    
    
    % % fast_sodium_current_j_gate
    if (V<-40.) 
         bbeta_j = ((0.02424*exp((-0.01052*V)))/(1.+exp((-0.1378*(V+40.14)))));
    else
       bbeta_j = ((0.6*exp((0.057*V)))/(1.+exp((-0.1*(V+32.))))); %units(unitless/ms);
    end
    j_inf = (1./((1.+exp(((V+71.55)/7.43)))*(1.+exp(((V+71.55)/7.43))))); %units(unitless);
    
    if (V<-40.)  
       aalpha_j = (((((-25428.*exp((0.2444*V))) - (6.948e-6*exp((-0.04391*V))))*(V+37.78))/1.)/(1.+exp((0.311*(V+79.23))))); 
    else
       aalpha_j =  0.; %units(unitless/ms);
    end
    tau_j = (1./(aalpha_j+bbeta_j)); %units(ms);
    diff_j = (j_inf - j)/tau_j;

% % sodium_background_current
    IbNa = (Gbna*(V - E_Na)); %units(pA/pF);

% % L_type_Ca_current
    ICal = ((((GCaL*d*f*f2*fCass*4.*(V - 15.)*(F*F))/(R*T))*((0.25*Cass*exp(((2.*(V - 15.)*F)/(R*T)))) - Cao))/(exp(((2.*(V - 15.)*F)/(R*T))) - 1.)); %units(pA/pF);

% %   L_type_Ca_current_d_gate
    aalpha_d = ((1.4/(1.+exp(((-35. - V)/13.))))+0.25); %units(unitless);
    gamma_d = (1./(1.+exp(((50. - V)/20.)))); %units(ms);
    bbeta_d = (1.4/(1.+exp(((V+5.)/5.)))); %units(unitless);
    tau_d = ((1.*aalpha_d*bbeta_d)+gamma_d); %units(ms);
    d_inf = (1./(1.+exp(((-8. - V)/7.5)))); %units(unitless);
    diff_d = (d_inf - d)/tau_d;

% %    L_type_Ca_current_f_gate
    tau_f = ((1102.5*exp((-((V+27.)*(V+27.))/225.)))+(200./(1.+exp(((13. - V)/10.))))+(180./(1.+exp(((V+30.)/10.))))+20.); %units(ms);
    f_inf = (1./(1.+exp(((V+20.)/7.)))); %units(unitless); 
    diff_f = ((f_inf - f)/tau_f); %units(unitless/ms);

% %  L_type_Ca_current_f2_gate
    tau_f2 = ((562.*exp((-((V+27.)*(V+27.))/240.)))+(31./(1.+exp(((25. - V)/10.))))+(80./(1.+exp(((V+30.)/10.))))); %units(ms);
    f2_inf = ((0.67/(1.+exp(((V+35.)/7.))))+0.33); %units(unitless);
    diff_f2 = ((f2_inf - f2)/tau_f2); %units(unitless/ms);   
       
% %   L_type_Ca_current_fCass_gate   
    fCass_inf = ((0.6/(1.+((Cass/0.05)*(Cass/0.05))))+0.4); %units(unitless);
    tau_fCass = ((80./(1.+((Cass/0.05)*(Cass/0.05))))+2.); %units(ms); 
    diff_fCass = ((fCass_inf - fCass)/tau_fCass); %units(unitless/ms);
  
% %  calcium_background_current
    IbCa = (Gbca*(V - E_Ca)); %units(pA/pF);  
    
% %    transient_outward_current
    Ito = (0.2*Gto*r*s*(V - E_K)); %units(pA/pF);    

% %     transient_outward_current_s_gate
    tau_s = ((85.*exp((-((V+25.)*(V+25.))/320.)))+(5./(1.+exp(((V - 40.)/5.))))+42.); %units(ms);
    s_inf = (1./(1.+exp(((V+27.)/13.)))); %units(unitless);
    diff_s = ((s_inf - s)/tau_s); %units(unitless/ms);

% % transient_outward_current_r_gate
    r_inf = (1./(1.+exp(((20. - V)/13.)))); %units(unitless);
    tau_r = ((10.45*exp((-((V+40.)*(V+40.))/1800.)))+7.3); %units(ms)
    diff_r = ((r_inf - r)/tau_r); %units(unitless/ms);

% %  sustained_outward_current
    a = (1./(1.+exp(((5. - V)/17.)))); %units(unitless);
    Isus = (Gsus*a*(V - E_K)); %units(pA/pF);   

% %   sodium_potassium_pump_current
    INaK = (((((maxINaK*Ko)/(Ko+K_mk))*Nai)/(Nai+K_mNa))/(1.+(0.1245*exp(((-0.1*V*F)/(R*T))))+(0.0353*exp(((-V*F)/(R*T)))))); %units(pA/pF);

% % sodium_calcium_exchanger_current
    INaCa = ((maxINaCa*((exp(((gamma*V*F)/(R*T)))*(Nai*Nai*Nai)*Cao) - (exp((((gamma - 1.)*V*F)/(R*T)))*(Nao*Nao*Nao)*Cai*alpha)))/(((Km_Nai*Km_Nai*Km_Nai)+...
    (Nao*Nao*Nao))*(Km_Ca+Cao)*(1.+(K_sat*exp((((gamma - 1.)*V*F)/(R*T))))))); %units(pA/pF);

% % calcium_pump_current
    IpCa = ((GpCa*Cai)/(Cai+K_pCa)); %units(pA/pF);

% % potassium_pump_current
    IpK = ((GpK*(V - E_K))/(1.+exp(((25. - V)/5.98)))); %units(pA/pF);

% %  calcium_dynamics
    Ca_sr_bufsr = (1./(1.+((Buf_sr*K_buf_sr)/((Ca_SR+K_buf_sr)*(Ca_SR+K_buf_sr))))); %units(unitless);
    Cass_bufss = (1./(1.+((Buf_ss*K_buf_ss)/((Cass+K_buf_ss)*(Cass+K_buf_ss))))); %units(unitless);
    kcasr = (max_sr - ((max_sr - min_sr)/(1.+((EC/Ca_SR)*(EC/Ca_SR))))); %units(unitless);
    k1 = (k1_prime/kcasr); %units(L^2/mmol^2/ms);
    O = ((k1*(Cass*Cass)*R_prime)/(k3+(k1*(Cass*Cass)))); %units(unitless);
    i_rel = (V_rel*O*(Ca_SR - Cass)); %units(mmol/L/ms);
    i_xfer = (V_xfer*(Cass - Cai)); %units(mmol/L/ms);

% %   calcium_dynamics Derivatives
    diff_Cass = (Cass_bufss*((((-1.*ICal*Cm)/(2.*1.*V_ss*F))+((i_rel*V_sr)/V_ss)) - ((i_xfer*V_c)/V_ss))); %units(mmol/L/ms);
    i_leak = (V_leak*(Ca_SR - Cai)); %units(mmol/L/ms);[tmat,sol]=
    i_up = (Vmax_up/(1.+((K_up*K_up)/(Cai*Cai)))); %units(mmol/L/ms);
    diff_Ca_SR = (Ca_sr_bufsr*(i_up - (i_rel+i_leak))); %units(mmol/L/ms);
    Cai_bufc = (1./(1.+((Buf_c*K_buf_c)/((Cai+K_buf_c)*(Cai+K_buf_c))))); %units(unitless);
    diff_Cai = (Cai_bufc*(((((i_leak - i_up)*V_sr)/V_c)+i_xfer) - ((1.*((IbCa+IpCa) - (2.*INaCa))*Cm)/(2.*1.*V_c*F)))); %units(mmol/L/ms);
    k2 = (k2_prime*kcasr); %units(L/mmol/ms);
    diff_R_prime = ((-k2*Cass*R_prime)+(k4*(1. - R_prime))); %units(unitless/ms);

% %     sodium_dynamics
    diff_Nai = (((-1.*(INa+IbNa+If_Na+(3.*INaK)+(3.*INaCa)))/(1.*V_c*F))*Cm); %units(mmol/L/ms);
    
% %     potassium_dynamics
    diff_Ki = (((-1.*((IK1+Ito+If_K+Isus+IKr+IKs+IpK) - (2.*INaK)))/(1.*V_c*F))*Cm); %units(mmol/L/ms);
    
% %  Stimulation   
   Istim = 0;
    for count=1:length(stimtimes)
        if (t >= stimtimes(count) && t < stimtimes(count)+stimDur)
            Istim = -60;
            break;
        end
        
    end


% %     Membrane current
    Iion = -(((-1./1.)*(IK1+Ito+Isus+IKr+IKs+ICal+INaK+INa+IbNa+INaCa+IbCa+IpK+IpCa+If))); %units(uA/uF);
    Itotal = (1/Cm)*(Istim +  Iion); 
    diff_V = -Itotal;
%%=============================================================================================================================
% % Purkinje Part of the simulation
%%==============================================================================================================================
% %     reversal_potentials
    E_Km = (((R*T)/F)*log((Ko/Kim))); %units(mV);
    E_Cam = (((0.5*R*T)/F)*log((Cao/Caim))); %units(mV);
    E_Ksm = (((R*T)/F)*log(((Ko+(P_kna*Nao))/(Kim+(P_kna*Naim))))); %units(mV);
    E_Nam = (((R*T)/F)*log((Nao/Naim))); %units(mV);

    
% %     h_inperpolarization_activated_current
    If_Nam =  ifMulp*(y_inm*Gf_Na*(Vm - E_Nam)); %units(pA/pF);
    If_Km =  ifMulp*(y_inm*Gf_K*(Vm - E_Km)); %units(pA/pF);
    Ifm =  ifMulp*(If_Nam+If_Km); %units(pA/pF);

% %     hyperpolarization_activated_current_y_gate
    y_infm = (1./(1.+exp(((Vm+80.6)/6.8)))); %units(unitless);
    aalpha_ym = (1.*exp((-2.9 - (0.04*Vm)))); %units(unitless/ms);
    bbeta_ym = (1.*exp((3.6+(0.11*Vm)))); %units(unitless/ms);
    tau_ym = (4000./(aalpha_ym+bbeta_ym)); %units(ms);
    diff_ym = (y_infm - y_inm)/tau_ym;

% % inward_rectifier_potassium_current
    xK1_infm = (1./(1.+exp((0.1*(Vm+75.44))))); %units(unitless);
    IK1m = (GK1*xK1_infm*((Vm - 8.) - E_Km)); %units(pA/pF);
    
% % rapid_time_dependent_potassium_current
    IKrm = (0.5*GKr*sqrt((Ko/5.4))*Xr1m*Xr2m*(Vm - E_Km)); %units(pA/pF);


% %     rapid_time_dependent_potassium_current_Xr1_gate
    xr1_infm = (1./(1.+exp(((-26. - Vm)/7.)))); %units(unitless);
    alpha_xr1m = (450./(1.+exp(((-45. - Vm)/10.)))); %units(unitless);
    beta_xr1m = (6./(1.+exp(((Vm+30.)/11.5)))); %units(unitless);
    tau_xr1m = (1.*alpha_xr1m*beta_xr1m); %units(ms);
    diff_Xr1m = ((xr1_infm - Xr1m)/tau_xr1m); %units(unitless/ms);

% %    rapid_time_dependent_potassium_current_Xr2_gate
    alpha_xr2m = (3./(1.+exp(((-60. - Vm)/20.)))); %units(unitless);
    xr2_infm = (1./(1.+exp(((Vm+88.)/24.)))); %units(unitless);
    beta_xr2m = (1.12/(1.+exp(((Vm - 60.)/20.)))); %units(unitless); 
    tau_xr2m = (1.*alpha_xr2m*beta_xr2m); %units(ms);
    diff_Xr2m = ((xr2_infm - Xr2m)/tau_xr2m); %units(unitless/ms);
    
% %     slow_time_dependent_potassium_current
    IKsm = (0.5*GKs*(Xsm*Xsm)*(Vm - E_Ksm)); %units(pA/pF);

% %     slow_time_dependent_potassium_current_Xs_gate
    alpha_xsm = (1400./sqrt((1.+exp(((5. - Vm)/6.))))); %units(unitless);
    beta_xsm = (1./(1.+exp(((Vm - 35.)/15.)))); %units(unitless);
    tau_xsm = ((1.*alpha_xsm*beta_xsm)+80.); %units(ms);
    xs_infm = (1./(1.+exp(((-5. - Vm)/14.)))); %units(unitless);
    diff_Xsm = ((xs_infm - Xsm)/tau_xsm); %units(unitless/ms);

% %     fast_sodium_current
    INam = (3*GNa*(mm*mm*mm)*hm*jm*(Vm - E_Nam)); %units(pA/pF);

% % fast_sodium_current_m_gate
    bbeta_mm = ((0.1/(1.+exp(((Vm+35.)/5.))))+(0.1/(1.+exp(((Vm - 50.)/200.))))); %units(unitless);
    aalpha_mm = (1./(1.+exp(((-60. - Vm)/5.)))); %units(unitless);
    m_infm = (1./((1.+exp(((-56.86 - Vm)/9.03)))*(1.+exp(((-56.86 - Vm)/9.03))))); %units(unitless);
    tau_mm = (1.*aalpha_mm*bbeta_mm); %units(ms);
    diff_mm = (m_infm - mm)/tau_mm;

% %     fast_sodium_current_h_gate
     if (Vm<-40.) 
       bbeta_hm = ((2.7*exp((0.079*Vm)))+(310000.*exp((0.3485*Vm)))); 
     else
          bbeta_hm =(0.77/(0.13*(1.+exp(((Vm+10.66)/-11.1))))); %units(unitless/ms);
     end
     if (Vm<-40.)
        aalpha_hm = (0.057*exp((-(Vm+80.)/6.8)));
     else
        aalpha_hm =  0.; %units(unitless/ms;
     end
    tau_hm = (1./(aalpha_hm+bbeta_hm)); %units(ms);
    h_infm = (1./((1.+exp(((Vm+71.55)/7.43)))*(1.+exp(((Vm+71.55)/7.43))))); %units(unitless);
    diff_hm = (h_infm - hm)/tau_hm;
    
    
    
    % % fast_sodium_current_j_gate
    if (Vm<-40.) 
         bbeta_jm = ((0.02424*exp((-0.01052*Vm)))/(1.+exp((-0.1378*(Vm+40.14)))));
    else
       bbeta_jm = ((0.6*exp((0.057*Vm)))/(1.+exp((-0.1*(Vm+32.))))); %units(unitless/ms);
    end
    j_infm = (1./((1.+exp(((Vm+71.55)/7.43)))*(1.+exp(((Vm+71.55)/7.43))))); %units(unitless);
    
    if (Vm<-40.)  
       aalpha_jm = (((((-25428.*exp((0.2444*Vm))) - (6.948e-6*exp((-0.04391*Vm))))*(Vm+37.78))/1.)/(1.+exp((0.311*(Vm+79.23))))); 
    else
       aalpha_jm =  0.; %units(unitless/ms);
    end
    tau_jm = (1./(aalpha_jm+bbeta_jm)); %units(ms);
    diff_jm = (j_infm - jm)/tau_jm;

% % sodium_background_current
    IbNam = (Gbna*(Vm - E_Nam)); %units(pA/pF);

% % L_type_Ca_current
    ICalm = ((((GCaL*dm*fm*f2m*fCassm*4.*(Vm - 15.)*(F*F))/(R*T))*((0.25*Cassm*exp(((2.*(Vm - 15.)*F)/(R*T)))) - Cao))/(exp(((2.*(Vm - 15.)*F)/(R*T))) - 1.)); %units(pA/pF);
    

% %   L_type_Ca_current_d_gate
    aalpha_dm = ((1.4/(1.+exp(((-35. - Vm)/13.))))+0.25); %units(unitless);
    gamma_dm = (1./(1.+exp(((50. - Vm)/20.)))); %units(ms);
    bbeta_dm = (1.4/(1.+exp(((Vm+5.)/5.)))); %units(unitless);
    tau_dm = ((1.*aalpha_dm*bbeta_dm)+gamma_dm); %units(ms);
    d_infm = (1./(1.+exp(((-8. - Vm)/7.5)))); %units(unitless);
    diff_dm = (d_infm - dm)/tau_dm;

% %    L_type_Ca_current_f_gate
    tau_fm = ((1102.5*exp((-((Vm+27.)*(Vm+27.))/225.)))+(200./(1.+exp(((13. - Vm)/10.))))+(180./(1.+exp(((Vm+30.)/10.))))+20.); %units(ms);
    f_infm = (1./(1.+exp(((Vm+20.)/7.)))); %units(unitless); 
    diff_fm = ((f_infm - fm)/tau_fm); %units(unitless/ms);

% %  L_type_Ca_current_f2_gate
    tau_f2m = ((562.*exp((-((Vm+27.)*(Vm+27.))/240.)))+(31./(1.+exp(((25. - Vm)/10.))))+(80./(1.+exp(((Vm+30.)/10.))))); %units(ms);
    f2_infm = ((0.67/(1.+exp(((Vm+35.)/7.))))+0.33); %units(unitless);
    diff_f2m = ((f2_infm - f2m)/tau_f2m); %units(unitless/ms);   
       
% %   L_type_Ca_current_fCass_gate   
    fCass_infm = ((0.6/(1.+((Cassm/0.05)*(Cassm/0.05))))+0.4); %units(unitless);
    tau_fCassm = ((80./(1.+((Cassm/0.05)*(Cassm/0.05))))+2.); %units(ms); 
    diff_fCassm = ((fCass_infm - fCassm)/tau_fCassm); %units(unitless/ms);
  
% %  calcium_background_current
    IbCam = (Gbca*(Vm - E_Cam)); %units(pA/pF);  
    
% %    transient_outward_current
    Itom = (2*Gto*rm*sm*(Vm - E_Km)); %units(pA/pF);    

% %     transient_outward_current_s_gate
    tau_sm = ((85.*exp((-((Vm+25.)*(Vm+25.))/320.)))+(5./(1.+exp(((Vm - 40.)/5.))))+42.); %units(ms);
    s_infm = (1./(1.+exp(((Vm+27.)/13.)))); %units(unitless);
    diff_sm = ((s_infm - sm)/tau_sm); %units(unitless/ms);

% % transient_outward_current_r_gate
    r_infm = (1./(1.+exp(((20. - Vm)/13.)))); %units(unitless);
    tau_rm = ((10.45*exp((-((Vm+40.)*(Vm+40.))/1800.)))+7.3); %units(ms)
    diff_rm = ((r_infm - rm)/tau_rm); %units(unitless/ms);

% %  sustained_outward_current
    am = (1./(1.+exp(((5. - Vm)/17.)))); %units(unitless);
    Isusm = (Gsus*am*(Vm - E_Km)); %units(pA/pF);   

% %   sodium_potassium_pump_current
    INaKm = (((((maxINaK*Ko)/(Ko+K_mk))*Naim)/(Naim+K_mNa))/(1.+(0.1245*exp(((-0.1*Vm*F)/(R*T))))+(0.0353*exp(((-Vm*F)/(R*T)))))); %units(pA/pF);

% % sodium_calcium_exchanger_current
    INaCam = ((maxINaCa*((exp(((gamma*Vm*F)/(R*T)))*(Naim*Naim*Naim)*Cao) - (exp((((gamma - 1.)*Vm*F)/(R*T)))*(Nao*Nao*Nao)*Caim*alpha)))/(((Km_Nai*Km_Nai*Km_Nai)+...
    (Nao*Nao*Nao))*(Km_Ca+Cao)*(1.+(K_sat*exp((((gamma - 1.)*Vm*F)/(R*T))))))); %units(pA/pF);

% % calcium_pump_current
    IpCam = ((GpCa*Caim)/(Caim+K_pCa)); %units(pA/pF);

% % potassium_pump_current
    IpKm = ((GpK*(Vm - E_Km))/(1.+exp(((25. - Vm)/5.98)))); %units(pA/pF);

% %  calcium_dynamics
    Ca_sr_bufsrm = (1./(1.+((Buf_sr*K_buf_sr)/((Ca_SRm+K_buf_sr)*(Ca_SRm+K_buf_sr))))); %units(unitless);
    Cass_bufssm = (1./(1.+((Buf_ss*K_buf_ss)/((Cassm+K_buf_ss)*(Cassm+K_buf_ss))))); %units(unitless);
    kcasrm = (max_sr - ((max_sr - min_sr)/(1.+((EC/Ca_SRm)*(EC/Ca_SRm))))); %units(unitless);
    k1m = (k1_prime/kcasrm); %units(L^2/mmol^2/ms);
    Om = ((k1m*(Cassm*Cassm)*R_primem)/(k3+(k1m*(Cassm*Cassm)))); %units(unitless);
    i_relm = (V_rel*Om*(Ca_SRm - Cassm)); %units(mmol/L/ms);
    i_xferm = (V_xfer*(Cassm - Caim)); %units(mmol/L/ms);

% %   calcium_dynamics Derivatives
    diff_Cassm = (Cass_bufssm*((((-1.*ICalm*Cm)/(2.*1.*V_ss*F))+((i_relm*V_sr)/V_ss)) - ((i_xferm*V_c)/V_ss))); %units(mmol/L/ms);
    i_leakm = (V_leak*(Ca_SRm - Caim)); %units(mmol/L/ms);[tmat,sol]=
    i_upm = (Vmax_up/(1.+((K_up*K_up)/(Caim*Caim)))); %units(mmol/L/ms);
    diff_Ca_SRm = (Ca_sr_bufsrm*(i_upm - (i_relm+i_leakm))); %units(mmol/L/ms);
    Cai_bufcm = (1./(1.+((Buf_c*K_buf_c)/((Caim+K_buf_c)*(Caim+K_buf_c))))); %units(unitless);
    diff_Caim = (Cai_bufcm*(((((i_leakm - i_upm)*V_sr)/V_c)+i_xferm) - ((1.*((IbCam+IpCam) - (2.*INaCam))*Cm)/(2.*1.*V_c*F)))); %units(mmol/L/ms);
    k2m = (k2_prime*kcasrm); %units(L/mmol/ms);
    diff_R_primem = ((-k2m*Cassm*R_primem)+(k4*(1. - R_primem))); %units(unitless/ms);

% %     sodium_dynamics
    diff_Naim = (((-1.*(INam+IbNam+If_Nam+(3.*INaKm)+(3.*INaCam)))/(1.*V_c*F))*Cm); %units(mmol/L/ms);
    
% %     potassium_dynamics
    diff_Kim = (((-1.*((IK1m+Itom+If_Km+Isusm+IKrm+IKsm+IpKm) - (2.*INaKm)))/(1.*V_c*F))*Cm); %units(mmol/L/ms);
    
% %  Stimulation   
   Istimp = 0;
% %     for count=1:length(stimtimes)
% %         if (t > stimtimes(count) && t < stimtimes(count)+stimDur)
% %             Istimp = -60;
% %         end
% %     end


% % Membrane current From ventricle
    Iionm = -(((-1./1.)*(IK1m+Itom+Isusm+IKrm+IKsm+ICalm+INaKm+INam+IbNam+INaCam+IbCam+IpKm+IpCam+Ifm))); %units(uA/uF);
    Itotalm = Istimp +  Iionm;
% %     diff_Vm = -Itotalm;
    diff_Vm = (1/Cm)*((1/K)*((V - Vm)/Rj) - Itotalm);
% %==============================================================================================================================
   
    out(1) = diff_V;
    out(2) = diff_y;
    out(3) = diff_Xr1;
    out(4) = diff_Xr2;
    out(5) = diff_Xs;
    out(6) = diff_m;
    out(7) = diff_h;
    out(8) = diff_j;
    out(9) = diff_d;
    out(10) = diff_f;
    out(11) = diff_f2;
    out(12) = diff_fCass;
    out(13) = diff_s;
    out(14) = diff_r;
    out(15) = diff_Cai;
    out(16) = diff_Ca_SR;
    out(17) = diff_Cass;
    out(18) = diff_R_prime;
    out(19) = diff_Nai;
    out(20) = diff_Ki;
% %     Start of Purkinje State Variables
    out(21) = diff_Vm;
    out(22) = diff_ym;
    out(23) = diff_Xr1m;
    out(24) = diff_Xr2m;
    out(25) = diff_Xsm;
    out(26) = diff_mm;
    out(27) = diff_hm;
    out(28) = diff_jm;
    out(29) = diff_dm;
    out(30) = diff_fm;
    out(31) = diff_f2m;
    out(32) = diff_fCassm;
    out(33) = diff_sm;
    out(34) = diff_rm;
    out(35) = diff_Caim;
    out(36) = diff_Ca_SRm;
    out(37) = diff_Cassm;
    out(38) = diff_R_primem;
    out(39) = diff_Naim;
    out(40) = diff_Kim;
    out = out';    
    
end
