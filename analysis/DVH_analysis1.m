clc;
clear;
close all;

%% set matRad runtime configuration
matRad_rc
param.logLevel=1;
%defaultRootPath = matRad_cfg.matRadRoot;
defaultRootPath = '\\compute-0-0\workspace';
job_folder='job1';
radiationMode='photons';
description='prostate';
caseID='3482'; % 3482 3648 3782 3790 3840
robustness_approach = 'nominal';
robustness='none'; % none COWC COWC2 c-COWC c-COWC2 INTERVAL2 INTERVAL3
plan_target='CTV'; % CTV PTV
plan_beams='9F';
plan_objectives='4';
shiftSD='5x10x5';
scen_mode='nomScen'; % nomScen impScen5 impScen_permuted_truncated5 impScen7 impScen_permuted_truncated7
wcFactor=1.0;
beta1=1/13;
beta2=13/13;
theta1=1.0;
theta2=0.1;

output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
    filesep plan_target filesep plan_beams filesep plan_objectives filesep shiftSD filesep scen_mode ];

%output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
%    filesep plan_target filesep plan_beams filesep plan_objectives filesep scen_mode filesep num2str(wcFactor) filesep num2str(beta1) '_to_' num2str(beta2) ];

%output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
%    filesep plan_target filesep plan_beams filesep plan_objectives filesep shiftSD filesep scen_mode filesep num2str(wcFactor) filesep num2str(beta1) '_to_' num2str(beta2) ];

%foldername = [defaultRootPath filesep '../../JOBS/cminimax2/1/job4' filesep output_folder];
foldername = [defaultRootPath filesep 'JOBS\cminimax2\artemisa\2023-03-23\1' filesep job_folder filesep output_folder];
listing = dir(foldername);
filename1=[foldername filesep listing(end).name filesep 'dvh_trustband_' robustness_approach '.fig'];
filename2=[foldername filesep listing(end).name filesep 'dvh_' robustness_approach '.fig'];
filename3=[foldername filesep listing(end).name filesep 'results.mat'];
filename4=[foldername filesep listing(end).name filesep 'plan.mat'];

%% Load results files
fig1=openfig(filename1);
fig2=openfig(filename2);
load(filename3);
load(filename4);

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

if isfield(run_config,'wcFactor')
    fprintf('wcFactor: \t %.2f \n', run_config.wcFactor);
end

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

%% Copying display names to patches
h1=findobj(fig1,'LineStyle',':');
h2=findobj(fig1,'FaceAlpha',0.2);
for structIx = 1:numel(h1)
    h2(structIx).DisplayName=h1(structIx).DisplayName;
end

%% Expected DVH
fprintf('!!!------------- Nominal DVH -------------!!!\n');

% Bladder
h=findobj(fig1,'LineStyle','-','Type','line','DisplayName','BLADDER');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x;y]);

f1=@(z) interp1(x,y,z);

x0=60; 
fprintf('BLADDER: \t V%.2f = %.2f %% \n', [x0;f1(x0)]);
dim = [0.4 0.4 .3 .3];
str = sprintf('V%.0f (Bladder) = %.2f %%', [x0;f1(x0)]);
annotation(fig2,'textbox',dim,'BackgroundColor','white','EdgeColor','black','FontSize',8,'String',str,'FitBoxToText','on');
%datatip(h(1),x0,f1(x0));

% Rectum
h=findobj(fig1,'LineStyle','-','Type','line','DisplayName','RECTUM');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x;y]);

f1=@(z) interp1(x,y,z);

x0=40; 
fprintf('RECTUM: \t V%.2f = %.2f %% \n', [x0;f1(x0)]);
dim = [0.3 0.5 .3 .3];
str = sprintf('V%.0f (Rectum) = %.2f %%', [x0;f1(x0)]);
annotation(fig2,'textbox',dim,'BackgroundColor','white','EdgeColor','black','FontSize',8,'String',str,'FitBoxToText','on');
%datatip(h(1),x0,f1(x0));

%% Expected DVH
fprintf('\n');
fprintf('!!!------------- c-DVH -------------!!!\n');

% Bladder
h=findobj(fig1,'LineStyle',':','Type','line','DisplayName','BLADDER');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x;y]);

f1=@(z) interp1(x,y,z);

h=findobj(fig1,'LineStyle','-','Type','patch','DisplayName','BLADDER');
x_min=h(1).XData(1:1000);
y_min=h(1).YData(1:1000);

f1_min=@(z) interp1(x_min,y_min,z);

x_max=h(1).XData(1601:3200);
y_max=h(1).YData(1601:3200);

f1_max=@(z) interp1(x_max,y_max,z);

x0=60;
fprintf('BLADDER: \t V%.2f = %.2f [ %.2f - %.2f ] %% \n', [x0;f1(x0);f1_min(x0);f1_max(x0)]);
dim = [0.5 0.1 .3 .3];
str = sprintf('V%.0f (Bladder) = %.2f [ %.2f - %.2f ] %%', [x0;f1(x0);f1_min(x0);f1_max(x0)]);
annotation(fig1,'textbox',dim,'BackgroundColor','white','EdgeColor','black','FontSize',8,'String',str,'FitBoxToText','on');
%datatip(h(1),x0,f1_max(x0));

% Rectum
h=findobj(fig1,'LineStyle',':','Type','line','DisplayName','RECTUM');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x;y]);

f1=@(z) interp1(x,y,z);

h=findobj(fig1,'LineStyle','-','Type','patch','DisplayName','RECTUM');
x_min=h(1).XData(1:1000);
y_min=h(1).YData(1:1000);

f1_min=@(z) interp1(x_min,y_min,z);

x_max=h(1).XData(1601:3200);
y_max=h(1).YData(1601:3200);

f1_max=@(z) interp1(x_max,y_max,z);

x0=40; 
fprintf('RECTUM: \t V%.2f = %.2f [ %.2f - %.2f ] %% \n', [x0;f1(x0);f1_min(x0);f1_max(x0)]);
dim = [0.4 0.2 .3 .3];
str = sprintf('V%.0f (Rectum) = %.2f [ %.2f - %.2f ] %%', [x0;f1(x0);f1_min(x0);f1_max(x0)]);
annotation(fig1,'textbox',dim,'BackgroundColor','white','EdgeColor','black','FontSize',8,'String',str,'FitBoxToText','on');
%datatip(h(1),x0,f1_max(x0));


%% CTV DVH
fprintf('\n');

% CTV
h=findobj(fig1,'LineStyle',':','Type','line','DisplayName','CTV');
x=h(1).XData;
y=h(1).YData;
%fprintf('%.2f, %.2f\n', [x1;y1]);

f1=@(z) interp1(x,y,z);

x0=78; 
fprintf('CTV: \t \t V%.2f = %.2f %% \n', [x0;f1(x0)]);


% CTV
fprintf('\t \t RI = %.4f \n', results.(['robustnessAnalysis_' robustness_approach]).robustnessIndex);
dim = [0.5 0.5 .3 .3];
str = sprintf('RI (CTV) = %.4f', results.(['robustnessAnalysis_' robustness_approach]).robustnessIndex);
annotation(fig1,'textbox',dim,'BackgroundColor','white','EdgeColor','black','FontSize',8,'String',str,'FitBoxToText','on');

% CTV
h=findobj(fig1,'LineStyle','-','Type','patch','DisplayName','CTV');
x1_tmp=h(1).XData(1:1000);
y1_tmp=h(1).YData(1:1000);

y1_first=find(y1_tmp==100,1,'last');
if(isempty(y1_first))
    y1_first=1;
end
y1_last=find(y1_tmp>=0,1,'last');
x1_tmp=x1_tmp(y1_first:y1_last);
y1_tmp=y1_tmp(y1_first:y1_last);
[y1,iy]=unique(y1_tmp);
x1=x1_tmp(iy);
f1=@(z) interp1(y1,x1,z, 'linear', 'extrap');

x3=linspace(0,100,100);
plot(f1(x3),x3,'black','LineStyle',':','LineWidth',2,'DisplayName','CTV Min');

x2_tmp=h(1).XData(1601:3200);
y2_tmp=h(1).YData(1601:3200);

y2_first=find(y2_tmp>=0,1,'first');
y2_last=find(y2_tmp==100,1,'first');
if(isempty(y2_last))
    y2_last=1;
end
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
delStructs={'PTV', 'RING 0 - 10 mm', 'RING 0 - 10 mm', 'RING 10 - 50 mm', 'BULB'};
for structIx=1:numel(delStructs)
    h1=findobj(fig1,'DisplayName',delStructs{structIx});
    for plotIx = 1:numel(h1)
        h1(plotIx).Visible=false;
    end
    h2=findobj(fig2,'DisplayName',delStructs{structIx});
    for plotIx = 1:numel(h2)
        h2(plotIx).Visible=false;
    end
end

fig1.Children(1).FontSize=7;
fig1.Position = [100 100 540 230];
fig2.Children(3).FontSize=7;

saveas(fig1,[foldername filesep listing(end).name filesep 'dvh_trustband_' robustness_approach],'png');
saveas(fig2,[foldername filesep listing(end).name filesep 'dvh_' robustness_approach],'png');

close all;