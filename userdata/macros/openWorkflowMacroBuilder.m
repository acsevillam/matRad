function openWorkflowMacroBuilder(varargin)
% openWorkflowMacroBuilder Open the planWorkflow GUI for macro creation.
% Usage: openWorkflowMacroBuilder('caseID','3482','rootPath',userDataRoot,'cacheRootPath',fullfile(userDataRoot,'output','cache'))
%
% This macro opens the interactive planWorkflow editor so templates and
% runnable workflow macros can be exported without executing workflow stages.

macroFolder = fileparts(mfilename('fullpath'));
macroRoot = fileparts(macroFolder);
if exist(fullfile(macroRoot,'helpers'),'dir') ~= 7
    macroRoot = macroFolder;
end
helperFolder = fullfile(macroRoot,'helpers');
if exist(helperFolder,'dir') == 7
    addpath(helperFolder);
end
if exist('ensurePlanWorkflowAvailable','file') == 2
    ensurePlanWorkflowAvailable('openWorkflowMacroBuilder');
end
userDataRoot = fileparts(macroRoot);

macroDefaults = struct();
macroDefaults.caseID = '3482';
macroDefaults.rootPath = userDataRoot;
macroDefaults.cacheRootPath = fullfile(userDataRoot,'output','cache');
macroOptions = planWorkflow.gui.PlanPresetWriter.parseMacroOptions(macroDefaults,varargin{:});

workflowConfig = struct();

% Shared configuration.
workflowConfig.rootPath = macroOptions.rootPath;
workflowConfig.outputRootPath = fullfile(macroOptions.rootPath,'output');
workflowConfig.patientDataPath = fullfile(macroOptions.rootPath,'patients');
workflowConfig.cacheRootPath = macroOptions.cacheRootPath;

% Initial editor defaults. These values are editable in the GUI.
workflowConfig.prepare.caseID = macroOptions.caseID;
workflowConfig.prepare.AcquisitionType = 'dicom';
workflowConfig.prepare.hlutFileName = 'matRad_default.hlut';
workflowConfig.prepare.description = 'prostate';
workflowConfig.prepare.plan_template = 'interval2_001';
workflowConfig.prepare.radiationMode = 'photons';
workflowConfig.prepare.plan_beams = '9F';
workflowConfig.prepare.dicomMetadata = struct();
workflowConfig.prepare.resolution = [3 3 3];
workflowConfig.prepare.n_cores = feature('numcores');

% Precompute stage.
workflowConfig.precompute.doseResolution = [3 3 3];
workflowConfig.precompute.reference.strategy = 'none';
workflowConfig.precompute.reference.scenario.mode = 'nomScen';
workflowConfig.precompute.reference.scenario.ctActive = false;
workflowConfig.precompute.reference.scenario.ctReferenceScenId = 1;

workflowConfig.precompute.robustPlans.robust_1.label = 'INTERVAL2';
workflowConfig.precompute.robustPlans.robust_1.objectiveSetName = 'robust_1';
workflowConfig.precompute.robustPlans.robust_1.strategy = 'INTERVAL2';
workflowConfig.precompute.robustPlans.robust_1.scenario.mode = 'wcScen';
workflowConfig.precompute.robustPlans.robust_1.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_1.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_1.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_1.scenario.wcSigma = 1.0;
workflowConfig.precompute.robustPlans.robust_1.variants.id = 'theta_10';
workflowConfig.precompute.robustPlans.robust_1.variants.label = 'theta1=10';
workflowConfig.precompute.robustPlans.robust_1.variants.theta1 = 10;

workflowConfig.precompute.useCache = true;
workflowConfig.precompute.writeCache = true;

% Dose-pulling stage.
workflowConfig.dosePulling.step1Enabled = false;
workflowConfig.dosePulling.step2Enabled = false;
workflowConfig.dosePulling.maxIterations = 100;
workflowConfig.dosePulling.scale_factor = 1;

% Optimize stage.
workflowConfig.optimize.optimizer = 'IPOPT';

% Sample stage.
workflowConfig.sampling.sampling_linkToOptimization = true;
workflowConfig.sampling.sampling_scen_mode = 'wcScen';
workflowConfig.sampling.sampling_ctActive = true;
workflowConfig.sampling.sampling_setupActive = true;
workflowConfig.sampling.sampling_rangeActive = false;
workflowConfig.sampling.sampling_shiftSD = [5 10 5];
workflowConfig.sampling.sampling_wcSigma = 1.5;

% Analyze stage.
workflowConfig.analysis.evaluationMode = 'total';
workflowConfig.analysis.doseWindow = [];
workflowConfig.analysis.doseWindowDvh = [];
workflowConfig.analysis.doseWindowUncertainty = [];
workflowConfig.analysis.doseWindowRelativeUncertainty1 = [];
workflowConfig.analysis.doseWindowRelativeUncertainty2 = [];
workflowConfig.analysis.doseWindowUvh = [];
workflowConfig.analysis.gammaWindow = [0 1];
workflowConfig.analysis.gammaCriteria = [3 3];
workflowConfig.analysis.robustnessCriteria = [5 5];
workflowConfig.analysis.robustnessTargetMode = 'include';
workflowConfig.analysis.robustnessTargets = {'CTV'};

workflow = planWorkflow.Workflow(workflowConfig);
cleanupObj = onCleanup(@() workflow.releaseMemory());

try
    workflow.gui();
catch ME
    if strcmp(ME.identifier,'planWorkflow:gui:PlanEditor:Cancelled')
        return;
    end
    rethrow(ME);
end

if ~isempty(workflow.guiProgressReporter)
    workflow.guiProgressReporter.log( ...
        'Macro builder finished. No workflow stages were executed.');
    workflow.guiProgressReporter.setProgress(0, ...
        'Macro builder finished.');
end

end
