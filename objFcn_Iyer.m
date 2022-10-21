function obj  = objFcn_Iyer(param)
% Function written by Anthony Owusu-Mensah
global target %% Target is the simulated IKr values from the original model
global stimtimes
% % global param
% % param(1) = 0.017147641733086;
% % param(2) = 0.03304608038835;
% % param(3) = 0.03969328381141;
% % param(4) = -0.04306054163980; 
% % param(5) = 0.02057448605977;
% % param(6) = 0.02617412715118;
% % param(7) = 0.00134366604423;
% % param(8) = -0.02691385498399;
% % param(9) = 0.10666316491288;
% % param(10) = 0.00568908859717;
% % param(11) = 0.00646393910049;
% % param(12) = -0.04536642959543;
% % param(13) = 0.00008039374403;
% % param(14) = 0.00000069808924;
% % param(15) = 0.1388; 
% % param(16) = 0.7891;

%Initial conditions, obtained after long term 1 Hz pacing
ic=[    -85.49190000, 
    0.00142380, 
    0.98646900,
    0.99907100,
    0.00142380,
    0.93722800,
    0.00036627,
    0.99816200,
    0.99816200,
    0.03065770,
    1.71221e-06,
    0.99764600,
    8.06835e-05,
    0.97941500,
    9.93905e-02,
    1.28293,
    0.000131111,
    0.0121366,
    0.0116031,
    0.132481,
    0.0981466,
    0.633156,
    0.00026248,
    0.04,
    2.67486e-07,
    1.94911e-06,
    0.879232,
    0.000336009,
    0.00244706,
    0.00243042,
    0.00996049,
    0.12297,
    0.136961,
    0.0079682,
    0.740619326616887,
    0.789728100000000,
    0.000679603250000000,
    0.000119228400000000,
    2.28398900000000e-05,
    5.80367300000000e-06


    ]; %VM

%time is measured in milliseconds.  May set duration of simulation here.
    tend=10000;

%By default, current clamp mode.  In vector stimtimes, you may set the
%times at which a stimulus current is applied.  The magnitude and duration
%of this stimulus is defined in fcn.m
    stimtimes=[1:1000:10000];

%If changes are made to the differential equations altering the stiffness
%of the problem, a maximum time step may be required for numerical
%integration.  
%ODE23T is selected given its performance characteristics,
%alternate options may be tried (help ode23T in the command window).
%On a modern PC, 0.25x real time simulation can be expected

%The solution is returned in sol; tmat is the times at which solutions are
%returned, the solutions_at determines the output times of solutions.
solutions_at=[0:1:tend];
% % tic;
opts=odeset('MaxStep', 0.5);
[tmat,sol]=ode23t(@(t,y)Aslanidi_Iyer(t,y,param), solutions_at, ic,opts);
OHerg = sol(:,39);
V = sol(:,1);
% % idx = find(tmat==5000);
% % toc;
T = 310.0;                       %units(K);
F     = 96486.7;        %units(J / kmol / K);
R     = 8314.3;         %units(C / mol);
RTonF = R * T / F;          %units(mV);
Ki      = 135.0;                 %units(mM);
Ke = 5.4;                        %units(mM);
E_K  = RTonF * log(Ke  / Ki );                                %units(mV);
GKr   = 0.015583333333333;       %units(mS / cm^2);
fKo = pow((Ke/4.0),0.5);
IKr = GKr.*fKo.*OHerg.*(V-E_K);
obj = sum((IKr-target).^2);
end


