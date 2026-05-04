function matRad_calcStudy(multScen,varargin)
% matRad uncertainty study wrapper
%
% call
%   matRad_calcStudy(multScen)
%   matRad_calcStudy(multScen,'SelectStructures',structSel)
%   matRad_calcStudy(multScen,'PatientMatFile',matPatientPath)
%
% input
%   multScen:           matRad scenario model used for uncertainty sampling
%   varargin:           optional name-value pairs
%                       'SelectStructures': structures which should be
%                       examined; empty selects all visible structures
%                       'OutputPath': folder for sampling output and report
%                       'PatientMatFile': absolute path to patient mat file;
%                       empty uses ct/cst/stf/pln/resultGUI from workspace
%                       'ListOfQI': quality indicator names for the report
%                       'OperatorName': operator name written to the report
%
% output
%   (binary)            all results are saved; a pdf report will be generated
%                       and saved
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

p = inputParser;
p.addRequired('multScen',@(x) isa(x,'matRad_ScenarioModel'));
p.addParameter('SelectStructures',cell(0),@iscellstr);
p.addParameter('OutputPath',matRad_cfg.userfolder{1},@isfolder);
p.addParameter('PatientMatFile','',@(x) isEmptyPath(x) || isfile(x));
p.addParameter('ListOfQI',{'mean', 'std', 'max', 'min', 'D_2', 'D_5', 'D_50', 'D_95', 'D_98'},@iscellstr);
p.addParameter('OperatorName','matRad User',@(x) isstring(x) || ischar(x));

%

p.parse(multScen,varargin{:});
multScen = p.Results.multScen;
outputPath = p.Results.OutputPath;
structSel = p.Results.SelectStructures;
matPatientPath = p.Results.PatientMatFile;
if isstring(matPatientPath) && isscalar(matPatientPath)
    matPatientPath = char(matPatientPath);
end
listOfQI = p.Results.ListOfQI;
operator = p.Results.OperatorName;


% require minimum number of scenarios to ensure proper statistics
if multScen.numScenarios() < 20
    matRad_cfg.dispWarning('Detected a low number of scenarios. Proceeding is not recommended.');
    sufficientStatistics = false;
    pause(1);
else
    sufficientStatistics = true;
end

%% load DICOM imported patient or run from workspace
if ~isempty(matPatientPath)
    patientData = load(matPatientPath,'ct','cst','stf','pln','resultGUI');
    requiredVars = {'ct','cst','stf','pln','resultGUI'};
    if ~all(isfield(patientData,requiredVars))
        matRad_cfg.dispError('PatientMatFile must contain ct, cst, stf, pln, and resultGUI.');
    end
    ct = patientData.ct;
    cst = patientData.cst;
    stf = patientData.stf;
    pln = patientData.pln;
    resultGUI = patientData.resultGUI;
else
    try
        ct          = evalin('base','ct');
        cst         = evalin('base','cst');
        stf         = evalin('base','stf');
        pln         = evalin('base','pln');
        resultGUI   = evalin('base','resultGUI');
    catch
        matRad_cfg.dispError('Workspace for sampling is incomplete.');
    end
end

% check if nominal workspace is complete
if ~(exist('ct','var') && exist('cst','var') && exist('stf','var') && exist('pln','var') && exist('resultGUI','var'))
    matRad_cfg.dispError('Workspace for sampling is incomplete.');
end

% calculate RBExDose
if ~isfield(pln, 'bioParam')
    if strcmp(pln.radiationMode, 'protons')
        pln.bioOptimization = 'RBExD';
        pln.model = 'constRBE';
    elseif strcmp(pln.radiationMode, 'carbon')
        pln.bioOptimization = 'RBExD';
        pln.model = 'LEM';
    end
    pln.bioParam = matRad_bioModel(pln.radiationMode, pln.bioOptimization, pln.model);
end


pln.robOpt   = false;
pln.sampling = true;

%% perform calculation and save
tic
[caSampRes, mSampDose, pln, resultGUInomScen]  = matRad_sampling(ct,stf,cst,pln,resultGUI.w,structSel,multScen);
computationTime = toc;

filename         = 'resultSampling';
save(filename, '-v7.3');

%% perform analysis
% start here loading resultSampling.mat if something went wrong during analysis or report generation
[structureStat, doseStat] = matRad_samplingAnalysis(ct,cst,pln,caSampRes,mSampDose,resultGUInomScen);

%% generate report
reportPath = 'report';
reportFolder = fullfile(outputPath,reportPath);
if ~exist(reportFolder,'dir')
    mkdir(reportFolder);
end

templatePath = fullfile(fileparts(mfilename('fullpath')),'main_template.tex');
if exist(templatePath,'file') ~= 2
    matRad_cfg.dispError('Could not find sampling report template: %s',templatePath);
end
copyfile(templatePath,fullfile(reportFolder,'main.tex'));
    
% generate actual latex report
success = matRad_latexReport(reportFolder,ct, cst, pln, resultGUInomScen, structureStat, doseStat, mSampDose, listOfQI,...
    'ComputationTime',computationTime,...
    'SufficientStatistics',sufficientStatistics,...
    'OperatorName',operator);


if success
    open(fullfile(reportFolder,'main.pdf'));
    
else
     matRad_cfg.dispError('Report PDF can not be opened...');
end

end

function tf = isEmptyPath(pathValue)

tf = isempty(pathValue) || (isstring(pathValue) && isscalar(pathValue) && strlength(pathValue) == 0);

end
