close all; clear; clc;
% % global param
param(1) = 0.017147641733086;
param(2) = 0.03304608038835;
param(3) = 0.03969328381141;
param(4) = -0.04306054163980; 
param(5) = 0.02057448605977;
param(6) = 0.02617412715118;
param(7) = 0.00134366604423;
param(8) = -0.02691385498399;
param(9) = 0.10666316491288;
param(10) = 0.00568908859717;
param(11) = 0.00646393910049;
param(12) = -0.04536642959543;
param(13) = 0.00008039374403;
param(14) = 0.00000069808924;
param(15) = 0.1388; 
param(16) = 0.7891;
param(17)	= 5.320000001;
param(18)   = 0.015583333333333;

global target
TraceFilename1 = '/home/cce2022/Desktop/AxelProject/ModifiedModels/Aslanid_org/benchSim/Trace_2.dat';
Asla = ReadTraceFile(TraceFilename1);
target = Asla(:,3);
x0 = param;
options = optimset('Display','iter','PlotFcns',@optimplotfval);
x = fminsearch(@obj_Iyer,x0,options); % Returned optimized parameters
% % nvars = 18;
% % x = ga(@objFcn_Iyer,nvars);
