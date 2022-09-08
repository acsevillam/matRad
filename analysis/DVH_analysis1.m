clc;
clear;
close all;

%% set matRad runtime configuration
matRad_rc
param.logLevel=1;
defaultRootPath = matRad_cfg.matRadRoot;
radiationMode='photons';
description='prostate';
caseID='3782'; % 3482 3648 3782
robustness='INTERVAL2'; % none COWC c-COWC
plan_target='CTV'; % CTV PTV
plan_beams='9F';
plan_objectives='4';
scen_mode='impScen_permuted_truncated'; % nomScen impScen impScen_permuted_truncated
wcFactor=1;
beta1=1/13;
beta2=13/13;
theta1=0.7;
%output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
%    filesep plan_target filesep plan_beams filesep plan_objectives filesep scen_mode filesep num2str(wcFactor) ]; %filesep '2022-07-10 03-59-36' filesep 'dvh_trustband_robust.fig'];

%output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
%    filesep plan_target filesep plan_beams filesep plan_objectives filesep scen_mode filesep num2str(wcFactor) filesep num2str(beta1) '_to_' num2str(beta2) ];

output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
    filesep plan_target filesep plan_beams filesep plan_objectives filesep scen_mode filesep num2str(wcFactor) filesep num2str(theta1) ];

foldername = [defaultRootPath filesep '../../JOBS/interval/1_all/job1' filesep output_folder];
listing = dir(foldername);
filename=[foldername filesep listing(end).name filesep 'dvh_trustband_robust.fig'];
openfig(filename);
%% Bladder
h=findobj(gca,'LineStyle',':','DisplayName','BLADDER');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x;y]);

f1=@(z) interp1(x,y,z);

x0=60; 
fprintf('BLADDER: \t V%.2f = %.2f\n', [x0;f1(x0)]);

%% Rectum
h=findobj(gca,'LineStyle',':','DisplayName','RECTUM');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x;y]);

f1=@(z) interp1(x,y,z);

x0=40; 
fprintf('RECTUM: \t V%.2f = %.2f\n', [x0;f1(x0)]);

%% CTV
h=findobj(gca,'LineStyle',':','DisplayName','CTV');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x1;y1]);

f1=@(z) interp1(x,y,z);

x0=78; 
fprintf('CTV: \t \t V%.2f = %.2f\n', [x0;f1(x0)]);

%% CTV
h=findobj(gca,'Type','line','DisplayName','CTV');
x1_tmp=h(2).XData;
y1_tmp=h(2).YData;

y1_first=find(y1_tmp==100,1,'last');
x1_tmp=x1_tmp(y1_first:end);
y1_tmp=y1_tmp(y1_first:end);
[y1,iy]=unique(y1_tmp);
x1=x1_tmp(iy);
f1=@(z) interp1(y1,x1,z);

x2_tmp=h(1).XData;
y2_tmp=h(1).YData;

%y2_first=find(y2_tmp==100,1,'last'); 
%x2_tmp=x2_tmp(y2_first+1:end);
%y2_tmp=y2_tmp(y2_first+1:end);
[y2,iy]=unique(y2_tmp);
x2=x2_tmp(iy);
f2=@(z) interp1(y2,x2,z);

x0=5; 
fprintf('CTV: \t \t V%.2f = [ %.2f - %.2f ]\n', [x0;f1(x0);f2(x0)]);

%%
%filename=[foldername filesep listing(end).name filesep 'dvh_robust.fig'];
%openfig(filename);