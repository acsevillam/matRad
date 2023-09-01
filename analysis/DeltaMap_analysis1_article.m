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
caseID='3790'; % 3482 3648 3782 3790 3840
robustness_approach = 'robust';
robustness='INTERVAL3'; % none COWC COWC2 c-COWC c-COWC2 INTERVAL2 INTERVAL3
plan_target='CTV'; % CTV PTV
plan_beams='9F';
plan_objectives='4';
shiftSD='5x10x5';
scen_mode='impScen5'; % nomScen impScen5 impScen_permuted_truncated5 impScen7 impScen_permuted_truncated7
wcFactor=1.0;
beta1=1/13;
p2=1;
beta2=p2/13;
theta1=1.0;
theta2=0.0;

%output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
%    filesep plan_target filesep plan_beams filesep plan_objectives filesep shiftSD filesep scen_mode ];

%output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
%    filesep plan_target filesep plan_beams filesep plan_objectives filesep shiftSD filesep scen_mode filesep num2str(wcFactor) filesep num2str(beta1) '_to_' num2str(beta2) ];

output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
    filesep plan_target filesep plan_beams filesep plan_objectives filesep shiftSD filesep scen_mode filesep num2str(wcFactor) filesep num2str(theta1) filesep num2str(theta2) ];


%foldername = [defaultRootPath filesep '../../JOBS/cminimax2/1/job4' filesep output_folder];
foldername = [defaultRootPath filesep 'JOBS\apolo\interval3\apolo\2023-08-30\1' filesep job_folder filesep output_folder];
listing = dir(foldername);
filename1=[foldername filesep listing(end).name filesep 'results.mat'];
filename2=[foldername filesep listing(end).name filesep 'plan.mat'];

%% Load results files
load(filename1);
load(filename2);

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
fprintf('Resolution: \t [%.2f %.2f %.2f] \n', [run_config.resolution(1),run_config.resolution(2),run_config.resolution(3)]);
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

%%

dotStructs={'PTV'};
for dotStructIx=1:numel(dotStructs)
    for  structIx = 1:size(cst,1)
        if(strcmp(cst{structIx,2},dotStructs{dotStructIx}))
            cst{structIx, 5}.visibleColor=[1,1,1];
        end
    end
end 

if strcmp(plan_target,'CTV')
    delStructs={'RING 0 - 20 mm','RING 20 - 50 mm','PTV'};
else
    delStructs={'RING 0 - 20 mm','RING 20 - 50 mm'};
end

for delStructIx=1:numel(delStructs)
    for  structIx = 1:size(cst,1)
        if(strcmp(cst{structIx,2},delStructs{delStructIx}))
            cst{structIx,3}='IGNORED';
        end
    end
end

% Create target mask
targetMask = zeros(size(results.(['robustnessAnalysis_' robustness_approach]).robustnessCube1));
for  i = 1:size(cst,1)
    if isequal(cst{i,3},'TARGET')
        targetMask(cst{i,4}{1,1}) = 1;
    end
end

ixBase=1;
ixRV=size(cst,1)+1;
cst{ixRV,1}=cst{end,1}+1;
cst{ixRV,2}='ROB DOSE';
cst{ixRV,3}='OAR';
cst{ixRV,4}{1,1}=find((results.(['robustnessAnalysis_' robustness_approach]).robustnessCube1<=1).*results.(['robustnessAnalysis_' robustness_approach]).robustnessCube1.*targetMask);
cst{ixRV,5}=cst{ixBase,5};
cst{ixRV,5}.visibleColor=[0,1,0];

ixBase=1;
ixRV=size(cst,1)+1;
cst{ixRV,1}=cst{end,1}+1;
cst{ixRV,2}='HIGH DOSE';
cst{ixRV,3}='OAR';
cst{ixRV,4}{1,1}=find((results.(['robustnessAnalysis_' robustness_approach]).robustnessCube1>=3).*results.(['robustnessAnalysis_' robustness_approach]).robustnessCube1.*targetMask);
cst{ixRV,5}=cst{ixBase,5};
cst{ixRV,5}.visibleColor=[0.75,0,0];

%isocenter = matRad_getIsoCenter(cst,ct,0);
%slice      = round(isocenter(3)./ct.resolution.z);
robustnessCriteria = [5 5];

for slice = 23:2:35
    [~, ~,robustnessFig1] = robustnessIndex1(results.(['robustnessAnalysis_' robustness_approach]).meanCubeW,results.(['robustnessAnalysis_' robustness_approach]).stdCubeW,results.(['robustnessAnalysis_' robustness_approach]).refDose,robustnessCriteria,slice,ct,cst);

    set(robustnessFig1,'PaperOrientation','landscape');
    set(robustnessFig1,'PaperPositionMode','auto');
    set(robustnessFig1.Children(4),'XLim',[76.5000 99.5000]);
    set(robustnessFig1.Children(4),'YLim',[70.5000 95.5000]);
    set(robustnessFig1.Children(4).Title,'String',strjoin(['interval3','theta1 =',string(theta1),', theta2 =',string(theta2)],' ')); % Nomimal PTV (Margin based) minimax strjoin(['cminimax',string(p2),'of 13'],' ') strjoin(['interval3','theta1:',string(theta1),' and theta2:',string(theta2)],' ')
    set(robustnessFig1.Children(6),'XLim',[76.5000 99.5000]);
    set(robustnessFig1.Children(6),'YLim',[70.5000 95.5000]);
    set(robustnessFig1.Children(6).Title,'String',strjoin(['interval3','theta1 =',string(theta1),', theta2 =',string(theta2)],' '));
    %set(robustnessFig1,'PaperSize',[5.8 4.0]);
    %print(robustnessFig1,[foldername filesep listing(end).name filesep 'robust_conformity2_' robustness_approach '_' plan_target '_' num2str(slice)],'-dpdf','-r0','-fillpage');
    %print(robustnessFig1,[foldername filesep listing(end).name filesep 'robust_conformity2_' robustness_approach '_' plan_target  '_' num2str(p2) '_13' '_' num2str(slice)],'-dpdf','-r0','-fillpage');
    print(robustnessFig1,[foldername filesep listing(end).name filesep 'robust_conformity2_' robustness_approach '_' plan_target  '_' num2str(theta1) '_' num2str(theta2) '_' num2str(slice)],'-dpdf','-r0','-fillpage');

end
close all;
diary off