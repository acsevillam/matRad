function runBreastMultipleWorkflowCore(prepareConfig,varargin)
% runBreastMultipleWorkflowCore Run the shared breast multiple-plan workflow.
% Usage: runBreastMultipleWorkflowCore(prepareConfig,'caseID','4136','rootPath',userDataRoot,'cacheRootPath',fullfile(userDataRoot,'output','cache'))

if nargin < 1 || ~isstruct(prepareConfig)
    error('planWorkflow:macros:InvalidPrepareConfig', ...
        'A scalar workflowConfig.prepare struct is required.');
end

macroFolder = fileparts(mfilename('fullpath'));
macroRoot = findMacroRoot(macroFolder);
helperFolder = fullfile(macroRoot,'helpers');
if exist(helperFolder,'dir') == 7
    addpath(helperFolder);
end
if exist('ensurePlanWorkflowAvailable','file') == 2
    ensurePlanWorkflowAvailable('runBreastMultipleWorkflowCore');
end

userDataRoot = fileparts(macroRoot);

macroDefaults = struct();
macroDefaults.caseID = fieldOrDefault(prepareConfig,'caseID','4136');
macroDefaults.rootPath = userDataRoot;
macroDefaults.cacheRootPath = fullfile(userDataRoot,'output','cache');
macroDefaults.randomSeed = [];
macroOptions = planWorkflow.gui.PlanPresetWriter.parseMacroOptions(macroDefaults,varargin{:});

workflowConfig = struct();

% Shared configuration.
workflowConfig.rootPath = macroOptions.rootPath;
workflowConfig.outputRootPath = fullfile(macroOptions.rootPath,'output');
workflowConfig.patientDataPath = fullfile(macroOptions.rootPath,'patients');
workflowConfig.cacheRootPath = macroOptions.cacheRootPath;

% Prepare stage.
workflowConfig.prepare = prepareConfig;
workflowConfig.prepare.caseID = macroOptions.caseID;

% Precompute stage.
workflowConfig.precompute.doseResolution = [5 5 5];
workflowConfig.precompute.reference.label = 'Nominal';
workflowConfig.precompute.reference.scenario.mode = 'nomScen';
workflowConfig.precompute.reference.scenario.ctActive = false;
workflowConfig.precompute.reference.scenario.ctReferenceScenId = 1;

workflowConfig.precompute.robustPlans.robust_1.label = 'PTV';
workflowConfig.precompute.robustPlans.robust_1.objectiveSetName = 'PTV';
workflowConfig.precompute.robustPlans.robust_1.scenario.mode = 'nomScen';
workflowConfig.precompute.robustPlans.robust_1.scenario.ctActive = false;
workflowConfig.precompute.robustPlans.robust_1.scenario.ctReferenceScenId = 1;

workflowConfig.precompute.robustPlans.robust_2.label = 'INTERVAL2';
workflowConfig.precompute.robustPlans.robust_2.objectiveSetName = 'Interval2';
workflowConfig.precompute.robustPlans.robust_2.dosePrecompute.useStreaming = true;
workflowConfig.precompute.robustPlans.robust_2.scenario.mode = 'impScen5';
workflowConfig.precompute.robustPlans.robust_2.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_2.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_2.scenario.shiftSD = [4 8 6];
workflowConfig.precompute.robustPlans.robust_2.scenario.wcSigma = 1.0;
workflowConfig.precompute.robustPlans.robust_2.variants(1).id = 'theta_10';
workflowConfig.precompute.robustPlans.robust_2.variants(1).label = 'theta1=10';
workflowConfig.precompute.robustPlans.robust_2.variants(1).theta1 = 10;

workflowConfig.precompute.useCache = true;
workflowConfig.precompute.writeCache = true;

% Dose-pulling stage.
workflowConfig.pullDose.step1Enabled = true;
workflowConfig.pullDose.step1Target = {'CTV'};
workflowConfig.pullDose.step1Criteria = {'COV1'};
workflowConfig.pullDose.step1Limit = 0.9;
workflowConfig.pullDose.step1Start = 10;
workflowConfig.pullDose.step2Enabled = false;
workflowConfig.pullDose.maxIterations = 100;

% Optimize stage.
workflowConfig.optimize.optimizer = 'IPOPT';

% Sample stage.
workflowConfig.sampling.sampling_linkToOptimization = true;
workflowConfig.sampling.sampling_scen_mode = 'impScen_permuted5_truncated';
workflowConfig.sampling.sampling_ctActive = true;
workflowConfig.sampling.sampling_setupActive = true;
workflowConfig.sampling.sampling_shiftSD = [4 8 6];
workflowConfig.sampling.sampling_wcSigma = 1.5;

% Analyze stage.
workflowConfig.analysis.evaluationMode = 'total';
workflowConfig.analysis.gammaWindow = [0 1];
workflowConfig.analysis.gammaCriteria = [3 3];
workflowConfig.analysis.robustnessCriteria = [5 5];
workflowConfig.analysis.robustnessTargetMode = 'include';
workflowConfig.analysis.robustnessTargets = {'CTV'};

workflow = planWorkflow.Workflow(workflowConfig);
cleanupObj = onCleanup(@() workflow.releaseMemory());

workflow.gui();
workflow.prepare();
workflow.precompute();
workflow.pullDose();
workflow.optimize();
workflow.sample();
workflow.analyze();

workflow.save();
workflow.releaseMemory();

end

function macroRoot = findMacroRoot(startFolder)

macroRoot = startFolder;
while true
    if exist(fullfile(macroRoot,'helpers'),'dir') == 7
        return;
    end
    parentFolder = fileparts(macroRoot);
    if strcmp(parentFolder,macroRoot)
        error('planWorkflow:macros:MacroRootNotFound', ...
            'Could not locate userdata/macros root from %s.',startFolder);
    end
    macroRoot = parentFolder;
end

end

function value = fieldOrDefault(s,fieldName,defaultValue)

if isfield(s,fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end

end
