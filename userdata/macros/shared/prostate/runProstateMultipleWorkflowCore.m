function runProstateMultipleWorkflowCore(prepareConfig,varargin)
% runProstateMultipleWorkflowCore Run the shared prostate multiple-plan workflow.
% Usage: runProstateMultipleWorkflowCore(prepareConfig,'caseID','3482','rootPath',userDataRoot,'cacheRootPath',fullfile(userDataRoot,'output','cache'),'randomSeed',[])

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
    ensurePlanWorkflowAvailable('runProstateMultipleWorkflowCore');
end
userDataRoot = fileparts(macroRoot);

macroDefaults = struct();
macroDefaults.caseID = fieldOrDefault(prepareConfig,'caseID','3482');
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
workflowConfig.precompute.doseResolution = [3 3 3];
workflowConfig.precompute.reference.label = 'Nominal';
workflowConfig.precompute.reference.scenario.mode = 'nomScen';
workflowConfig.precompute.reference.scenario.ctActive = false;
workflowConfig.precompute.reference.scenario.ctReferenceScenId = 1;

workflowConfig.precompute.robustPlans.robust_1.label = 'PTV';
workflowConfig.precompute.robustPlans.robust_1.objectiveSetName = 'PTV';
workflowConfig.precompute.robustPlans.robust_1.scenario.mode = 'nomScen';
workflowConfig.precompute.robustPlans.robust_1.scenario.ctActive = false;
workflowConfig.precompute.robustPlans.robust_1.scenario.ctReferenceScenId = 1;

workflowConfig.precompute.robustPlans.robust_2.label = 'Minimax';
workflowConfig.precompute.robustPlans.robust_2.objectiveSetName = 'Minimax';
workflowConfig.precompute.robustPlans.robust_2.scenario.mode = 'impScen5';
workflowConfig.precompute.robustPlans.robust_2.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_2.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_2.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_2.scenario.wcSigma = 1.0;

workflowConfig.precompute.robustPlans.robust_3.label = 'Stochastic';
workflowConfig.precompute.robustPlans.robust_3.objectiveSetName = 'Stochastic';
workflowConfig.precompute.robustPlans.robust_3.scenario.mode = 'impScen5';
workflowConfig.precompute.robustPlans.robust_3.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_3.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_3.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_3.scenario.wcSigma = 1.0;

workflowConfig.precompute.robustPlans.robust_4.label = 'c-Minimax';
workflowConfig.precompute.robustPlans.robust_4.objectiveSetName = 'cMinimax';
workflowConfig.precompute.robustPlans.robust_4.scenario.mode = 'impScen5';
workflowConfig.precompute.robustPlans.robust_4.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_4.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_4.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_4.scenario.wcSigma = 1.0;
workflowConfig.precompute.robustPlans.robust_4.variants(1).id = 'p1_1_p2_2';
workflowConfig.precompute.robustPlans.robust_4.variants(1).label = 'p1=1 - p2=2';
workflowConfig.precompute.robustPlans.robust_4.variants(1).p1 = 1;
workflowConfig.precompute.robustPlans.robust_4.variants(1).p2 = 2;

workflowConfig.precompute.robustPlans.robust_5.label = 'MeanVariance';
workflowConfig.precompute.robustPlans.robust_5.objectiveSetName = 'MeanVariance';
workflowConfig.precompute.robustPlans.robust_5.scenario.mode = 'impScen5';
workflowConfig.precompute.robustPlans.robust_5.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_5.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_5.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_5.scenario.wcSigma = 1.0;

workflowConfig.precompute.robustPlans.robust_10.label = 'INTERVAL2';
workflowConfig.precompute.robustPlans.robust_10.objectiveSetName = 'Interval2';
workflowConfig.precompute.robustPlans.robust_10.scenario.mode = 'impScen5';
workflowConfig.precompute.robustPlans.robust_10.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_10.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_10.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_10.scenario.wcSigma = 1.0;
workflowConfig.precompute.robustPlans.robust_10.variants(1).id = 'theta1_1';
workflowConfig.precompute.robustPlans.robust_10.variants(1).label = 'theta1=1';
workflowConfig.precompute.robustPlans.robust_10.variants(1).theta1 = 1;
workflowConfig.precompute.robustPlans.robust_10.variants(2).id = 'theta1_2';
workflowConfig.precompute.robustPlans.robust_10.variants(2).label = 'theta1=2';
workflowConfig.precompute.robustPlans.robust_10.variants(2).theta1 = 2;
workflowConfig.precompute.robustPlans.robust_10.variants(3).id = 'theta1_5';
workflowConfig.precompute.robustPlans.robust_10.variants(3).label = 'theta1=5';
workflowConfig.precompute.robustPlans.robust_10.variants(3).theta1 = 5;
workflowConfig.precompute.robustPlans.robust_10.variants(4).id = 'theta1_10';
workflowConfig.precompute.robustPlans.robust_10.variants(4).label = 'theta1=10';
workflowConfig.precompute.robustPlans.robust_10.variants(4).theta1 = 10;
workflowConfig.precompute.robustPlans.robust_10.variants(5).id = 'theta1_20';
workflowConfig.precompute.robustPlans.robust_10.variants(5).label = 'theta1=20';
workflowConfig.precompute.robustPlans.robust_10.variants(5).theta1 = 20;
workflowConfig.precompute.robustPlans.robust_10.variants(6).id = 'theta1_0p01';
workflowConfig.precompute.robustPlans.robust_10.variants(6).label = 'theta1=0.01';
workflowConfig.precompute.robustPlans.robust_10.variants(6).theta1 = 0.01;
workflowConfig.precompute.robustPlans.robust_10.variants(7).id = 'theta1_0p02';
workflowConfig.precompute.robustPlans.robust_10.variants(7).label = 'theta1=0.02';
workflowConfig.precompute.robustPlans.robust_10.variants(7).theta1 = 0.02;
workflowConfig.precompute.robustPlans.robust_10.variants(8).id = 'theta1_0p05';
workflowConfig.precompute.robustPlans.robust_10.variants(8).label = 'theta1=0.05';
workflowConfig.precompute.robustPlans.robust_10.variants(8).theta1 = 0.05;
workflowConfig.precompute.robustPlans.robust_10.variants(9).id = 'theta1_0p1';
workflowConfig.precompute.robustPlans.robust_10.variants(9).label = 'theta1=0.1';
workflowConfig.precompute.robustPlans.robust_10.variants(9).theta1 = 0.1;
workflowConfig.precompute.robustPlans.robust_10.variants(10).id = 'theta1_0p2';
workflowConfig.precompute.robustPlans.robust_10.variants(10).label = 'theta1=0.2';
workflowConfig.precompute.robustPlans.robust_10.variants(10).theta1 = 0.2;
workflowConfig.precompute.robustPlans.robust_10.variants(11).id = 'theta1_0p5';
workflowConfig.precompute.robustPlans.robust_10.variants(11).label = 'theta1=0.5';
workflowConfig.precompute.robustPlans.robust_10.variants(11).theta1 = 0.5;
workflowConfig.precompute.robustPlans.robust_10.variants(12).id = 'theta1_50';
workflowConfig.precompute.robustPlans.robust_10.variants(12).label = 'theta1=50';
workflowConfig.precompute.robustPlans.robust_10.variants(12).theta1 = 50;

workflowConfig.precompute.robustPlans.robust_11.label = 'INTERVAL3';
workflowConfig.precompute.robustPlans.robust_11.objectiveSetName = 'Interval3';
workflowConfig.precompute.robustPlans.robust_11.scenario.mode = 'impScen5';
workflowConfig.precompute.robustPlans.robust_11.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_11.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_11.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_11.scenario.wcSigma = 1.0;
workflowConfig.precompute.robustPlans.robust_11.variants(1).id = 'theta1_1_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(1).label = 'theta1=1 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(1).theta1 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(1).theta2 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(2).id = 'theta1_2_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(2).label = 'theta1=2 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(2).theta1 = 2;
workflowConfig.precompute.robustPlans.robust_11.variants(2).theta2 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(3).id = 'theta1_5_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(3).label = 'theta1=5 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(3).theta1 = 5;
workflowConfig.precompute.robustPlans.robust_11.variants(3).theta2 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(4).id = 'theta1_10_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(4).label = 'theta1=10 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(4).theta1 = 10;
workflowConfig.precompute.robustPlans.robust_11.variants(4).theta2 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(5).id = 'theta1_20_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(5).label = 'theta1=20 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(5).theta1 = 20;
workflowConfig.precompute.robustPlans.robust_11.variants(5).theta2 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(6).id = 'theta1_0p01_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(6).label = 'theta1=0.01 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(6).theta1 = 0.01;
workflowConfig.precompute.robustPlans.robust_11.variants(6).theta2 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(7).id = 'theta1_0p02_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(7).label = 'theta1=0.02 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(7).theta1 = 0.02;
workflowConfig.precompute.robustPlans.robust_11.variants(7).theta2 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(8).id = 'theta1_0p05_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(8).label = 'theta1=0.05 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(8).theta1 = 0.05;
workflowConfig.precompute.robustPlans.robust_11.variants(8).theta2 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(9).id = 'theta1_0p1_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(9).label = 'theta1=0.1 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(9).theta1 = 0.1;
workflowConfig.precompute.robustPlans.robust_11.variants(9).theta2 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(10).id = 'theta1_0p2_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(10).label = 'theta1=0.2 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(10).theta1 = 0.2;
workflowConfig.precompute.robustPlans.robust_11.variants(10).theta2 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(11).id = 'theta1_0p5_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(11).label = 'theta1=0.5 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(11).theta1 = 0.5;
workflowConfig.precompute.robustPlans.robust_11.variants(11).theta2 = 1;
workflowConfig.precompute.robustPlans.robust_11.variants(12).id = 'theta1_50_theta2_1';
workflowConfig.precompute.robustPlans.robust_11.variants(12).label = 'theta1=50 - theta2=1';
workflowConfig.precompute.robustPlans.robust_11.variants(12).theta1 = 50;
workflowConfig.precompute.robustPlans.robust_11.variants(12).theta2 = 1;

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
workflowConfig.sampling.sampling_shiftSD = [5 10 5];
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
