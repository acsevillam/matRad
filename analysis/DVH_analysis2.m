clc;
clear;
close all;

%% set matRad runtime configuration
matRad_rc
param.logLevel=1;
defaultRootPath = matRad_cfg.matRadRoot;
radiationMode='photons';
description='breast';
caseID='3832'; % 3477 3749 3832 3833 3929
robustness='INTERVAL2'; % none COWC c-COWC INTERVAL2 INTERVAL3
plan_target='CTV'; % CTV PTV
plan_beams='7F';
plan_objectives='4';
scen_mode='impScen_permuted_truncated'; % nomScen impScen impScen_permuted_truncated
wcFactor=1.0;
beta1=1/13;
beta2=13/13;
theta1=0.8;
theta2=0.1;
%output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
%    filesep plan_target filesep plan_beams filesep plan_objectives filesep scen_mode filesep num2str(wcFactor) ]; %filesep '2022-07-10 03-59-36' filesep 'dvh_trustband_robust.fig'];

%output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
%    filesep plan_target filesep plan_beams filesep plan_objectives filesep scen_mode filesep num2str(wcFactor) filesep num2str(beta1) '_to_' num2str(beta2) ];

output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
    filesep plan_target filesep plan_beams filesep plan_objectives filesep scen_mode filesep num2str(wcFactor) filesep num2str(theta1) ];

%output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
%    filesep plan_target filesep plan_beams filesep plan_objectives filesep scen_mode filesep num2str(wcFactor) filesep num2str(theta1) filesep num2str(theta2) ];

foldername = [defaultRootPath filesep '../../JOBS/interval2/2_all/job1' filesep output_folder];
listing = dir(foldername);
filename1=[foldername filesep listing(end).name filesep 'dvh_trustband_robust.fig'];
filename2=[foldername filesep listing(end).name filesep 'results.mat'];
filename3=[foldername filesep listing(end).name filesep 'plan.mat'];

%%
openfig(filename1);
load(filename2);
load(filename3);

%% Initiallize diary log
diaryname=[foldername filesep listing(end).name filesep 'results.log'];

if exist(diaryname, 'file')==2
  delete(diaryname);
end
diary(diaryname)
diary on

%% Description
fprintf('Description: \t %s \n', run_config.description);
fprintf('CaseID: \t %s \n', run_config.caseID);
fprintf('Resolution: \t %s \n', run_config.resolution);
fprintf('Objectives: \t %s \n', run_config.plan_objectives);
fprintf('Target: \t %s \n', run_config.plan_target);
fprintf('Beam setup: \t %s \n', run_config.plan_beams);
fprintf('shiftSD: \t [%.2f %.2f %.2f] \n', [run_config.shiftSD(1),run_config.shiftSD(2),run_config.shiftSD(3)]);
fprintf('Robustness: \t %s \n', run_config.robustness);
fprintf('Scen mode: \t %s \n', run_config.scen_mode);
fprintf('wcFactor: \t %.2f \n', run_config.wcFactor);

if isfield(run_config,'beta1')
    fprintf('beta1: \t %.2f \n', run_config.beta1);
end

if isfield(run_config,'beta2')
    fprintf('beta2: \t %.2f \n', run_config.beta2);
end

if isfield(run_config,'theta1')
    fprintf('theta1: \t %.2f \n', run_config.theta1);
end

if isfield(run_config,'theta2')
    fprintf('theta2: \t %.2f \n', run_config.theta2);
end

fprintf('Samp. mode: \t %s \n', run_config.sampling_mode);
fprintf('Samp. wcFactor:  %.2f \n\n\n', run_config.sampling_wcFactor);

%% Left Lung
h=findobj(gca,'LineStyle',':','DisplayName','LEFT LUNG');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x;y]);

f1=@(z) interp1(x,y,z);

x0=20; 
fprintf('LEFT LUNG: \t V%.2f = %.2f %% \n', [x0;f1(x0)]);

%% Heart
h=findobj(gca,'LineStyle',':','DisplayName','HEART');
x=h(1).XData;
y=h(1).YData;
step=x(2)-x(1);
d = -diff(y)/step;
x=x(1:end-1);
%fprintf('%.2f, %.2f\n', [x;y]);

f1=@(z) interp1(x,d,z);
mean_dose = x*f1(x)'/sum(f1(x));

fprintf('HEART: \t DMEAN = %.2f Gy \n', mean_dose);

%% CTV
h=findobj(gca,'LineStyle',':','DisplayName','CTV');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x1;y1]);

f1=@(z) interp1(x,y,z);

x0=42.56; 
fprintf('CTV: \t \t V%.2f = %.2f %% \n', [x0;f1(x0)]);

%% CTV
fprintf('\t \t RI = %.4f \n', results.robustnessAnalysis_robust.robustnessIndex);

%% CTV
h=findobj(gca,'LineStyle','-','FaceColor',[1,0,1]);
x1_tmp=h(1).XData(1:1000);
y1_tmp=h(1).YData(1:1000);

y1_first=find(y1_tmp==100,1,'last');
y1_last=find(y1_tmp>=0,1,'last');
x1_tmp=x1_tmp(y1_first:y1_last);
y1_tmp=y1_tmp(y1_first:y1_last);
[y1,iy]=unique(y1_tmp);
x1=x1_tmp(iy);
f1=@(z) interp1(y1,x1,z, 'linear', 'extrap');

x3=linspace(0,100,100);
plot(f1(x3),x3,'black','LineStyle',':','LineWidth',2,'DisplayName','CTV Min');

x2_tmp=h(1).XData(1001:2000);
y2_tmp=h(1).YData(1001:2000);

y2_first=find(y2_tmp>=0,1,'first');
y2_last=find(y2_tmp==100,1,'first');

%y2_first=find(y2_tmp==100,1,'last'); 
x2_tmp=x2_tmp(y2_first+1:y2_last);
y2_tmp=y2_tmp(y2_first+1:y2_last);
[y2,iy]=unique(y2_tmp);
x2=x2_tmp(iy);
f2=@(z) interp1(y2,x2,z, 'linear', 'extrap');

x0=5;
fprintf('\t \t Tail width (V%.2f) = %.2f [ %.2f - %.2f ] Gy \n', [x0;f2(x0)-f1(x0);f1(x0);f2(x0)]);

plot(f2(x3),x3,'black','LineStyle',':','LineWidth',2,'DisplayName','CTV Max');

%%
diary off

%%
%filename=[foldername filesep listing(end).name filesep 'dvh_robust.fig'];
%openfig(filename);
