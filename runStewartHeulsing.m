close all; clear; clc
%Script is written by Anthony Owusu-Mensah
% This executes function StewartHeulsing
% and generates the resulting action potential from
% Ventricular and Purkinje. 
%Initial conditions, obtained after long term 1 Hz pacing
ic=[  -69.1370441635924,
    0.0457562667986602,
    0.00550281999719088,
    0.313213286437995,
    0.00953708522974789,
    0.0417391656294997,
    0.190678733735145,
    0.238219836154029,
    0.000287906256206415,
    0.989328560287987,
    0.995474890442185,
    0.999955429598213,
    0.96386101799501,
    0.00103618091196912,
    0.000101878186157052,
    3.10836886659417,
    0.000446818714055411,
    0.991580051907845,
    8.80420286531673,
    136.781894160227,
    -69.1370441635924,
    0.0457562667986602,
    0.00550281999719088,
    0.313213286437995,
    0.00953708522974789,
    0.0417391656294997,
    0.190678733735145,
    0.238219836154029,
    0.000287906256206415,
    0.989328560287987,
    0.995474890442185,
    0.999955429598213,
    0.96386101799501,
    0.00103618091196912,
    0.000101878186157052,
    3.10836886659417,
    0.000446818714055411,
    0.991580051907845,
    8.80420286531673,
    136.781894160227
    ];
% % ic = ic';
tend=10000;
stimtimes= 0:1000:10000;
solutions_at=0:1:tend;
tic;
opts=odeset('MaxStep', 0.5);
[tmat,sol]= ode23t(@(t,y)StewartHeulsing(t,y,stimtimes), solutions_at, ic,opts);
toc;
figure('Color','w')
plot(tmat,sol(:,1) ,'b','LineWidth',1.5); hold on;
plot(tmat,sol(:,21) ,'r','LineWidth',1.5);
legend("VM","PC");
xlim([6900 tmat(end)])
title('Ventricle => Purkinje conduction')
% % plot(tmat,sol(:,1) ,'b','LineWidth',1.5);