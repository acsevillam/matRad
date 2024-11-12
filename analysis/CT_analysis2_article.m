clc;
clear;
close all;

%% set matRad runtime configuration
matRad_rc
param.logLevel=1;
%defaultRootPath = matRad_cfg.matRadRoot;
defaultRootPath = '/Volumes/BACKUP/workspace';
job_folder='job1';
radiationMode='photons';
description='breast';
caseID='4585'; % 3477 3749 3832 3833 3929 4136 4155 4203 4357 4390 4428 4494 4531 4585 4681
robustness_approach = 'nominal';
robustness='none'; % none COWC COWC2 c-COWC c-COWC2 INTERVAL2 INTERVAL3
plan_target='CTV'; % CTV PTV
plan_beams='7F';
plan_objectives='4';
shiftSD='4x8x6';
scen_mode='nomScen'; % nomScen impScen5 impScen_permuted_truncated5 impScen7 impScen_permuted_truncated7
wcFactor=1.0;
beta1=1/13;
p2=1;
beta2=p2/13;
theta1=1.0;
theta2=0.1;

output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
    filesep plan_target filesep plan_beams filesep plan_objectives filesep shiftSD filesep scen_mode ];

%output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
%    filesep plan_target filesep plan_beams filesep plan_objectives filesep scen_mode filesep num2str(wcFactor) filesep num2str(beta1) '_to_' num2str(beta2) ];

%output_folder = ['output' filesep radiationMode filesep description filesep caseID filesep robustness ...
%    filesep plan_target filesep plan_beams filesep plan_objectives filesep shiftSD filesep scen_mode filesep num2str(wcFactor) filesep num2str(beta1) '_to_' num2str(beta2) ];

%foldername = [defaultRootPath filesep '../../JOBS/cminimax2/1/job4' filesep output_folder];
%foldername = [defaultRootPath filesep 'JOBS\apolo\cminimax2\2023-06-11\2' filesep job_folder filesep output_folder];
foldername = [defaultRootPath filesep 'JOBS/apolo/cminimax2/2024-10-10/2' filesep job_folder filesep output_folder];
listing = dir(foldername);
filename1=[foldername filesep listing(end).name filesep 'ct.fig'];
filename2=[foldername filesep listing(end).name filesep 'dvh_' robustness_approach '.fig'];

%% Load results files
fig1=openfig(filename1);
fig2=openfig(filename2);

%%
lines = findobj(fig1,'Type','line');

for lIx=1:numel(lines)
  lines(lIx).LineWidth = 2;
end

%%
delStructs={'RING 0 - 20 mm', 'RING 0 - 20 mm', 'RING 20 - 50 mm', 'CTV Min', 'CTV Max'};
for structIx=1:numel(delStructs)
  l=findobj(fig2,'LineStyle','-','Type','line','DisplayName',delStructs{structIx});
  if exist('l') && ~isempty(l)
    delete(findobj(fig1,'Type','line','Color',l.Color));
  end
end

%%
colStructs={'PTV','SKIN'};
col={[1 0 0],[0 0 1]};
for structIx=1:numel(colStructs)
  l=findobj(fig2,'LineStyle','-','Type','line','DisplayName',colStructs{structIx});
  if exist('l') && ~isempty(l)
    lines=findobj(fig1,'Type','line','Color',l.Color);
    for lIx=1:numel(lines)
      lines(lIx).Color=col{structIx};
    end
  end
end

%%
fig1.Children(1).FontSize=9;
fig1.Position = [120 100 540 230];
fig1.Children(2).Title.String= 'Nomimal CTV (Margin based)';%strjoin(['cminimax',string(p2),'of 13'],' '); % Nomimal PTV (Margin based) minimax strjoin(['cminimax',string(p2),'of 13'],' ')

fig2.Children(4).Title.String= 'Nomimal CTV (Margin based)';%strjoin(['cminimax',string(p2),'of 13'],' '); % Nomimal PTV (Margin based) minimax strjoin(['cminimax',string(p2),'of 13'],' ')
fig2.Children(3).FontSize=9;

%saveas(fig1,[foldername filesep listing(end).name filesep 'dvh_trustband_' robustness_approach '_2'],'png');
set(fig1,'PaperOrientation','portrait');
set(fig1,'PaperPositionMode','auto');
set(fig1,'PaperSize',[4.0 8.0]);
%print(fig1,[foldername filesep listing(end).name filesep 'dvh_trustband_' robustness_approach '_' num2str(p2) '_13'],'-dpdf','-r0','-fillpage');
print(fig1,[foldername filesep listing(end).name filesep 'ct'],'-dpdf','-r0','-fillpage');

close all;