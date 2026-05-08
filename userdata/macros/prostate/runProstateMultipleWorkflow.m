function runProstateMultipleWorkflow(varargin)
% runProstateMultipleWorkflow Run the prostate multiple-plan workflow.
% Usage: runProstateMultipleWorkflow('caseID','3482','rootPath',userDataRoot,'cacheRootPath',fullfile(userDataRoot,'output','cache'),'randomSeed',[])

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
    ensurePlanWorkflowAvailable('runProstateMultipleWorkflow');
end
userDataRoot = fileparts(macroRoot);

macroDefaults = struct();
macroDefaults.caseID = '3482';
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
workflowConfig.prepare.caseID = macroOptions.caseID;
workflowConfig.prepare.AcquisitionType = 'dicom';
workflowConfig.prepare.hlutFileName = 'matRad_default.hlut';
workflowConfig.prepare.description = 'prostate';
workflowConfig.prepare.plan_template = 'comparison_001';
workflowConfig.prepare.radiationMode = 'photons';
workflowConfig.prepare.machine = 'Generic';
workflowConfig.prepare.bioModel = 'none';
workflowConfig.prepare.plan_beams = '9F';
workflowConfig.prepare.dicomMetadata = struct();
workflowConfig.prepare.resolution = [3 3 3];
workflowConfig.prepare.n_cores = 8;

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
workflowConfig.precompute.robustPlans.robust_2.scenario.mode = 'wcScen';
workflowConfig.precompute.robustPlans.robust_2.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_2.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_2.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_2.scenario.wcSigma = 1.0;

workflowConfig.precompute.robustPlans.robust_3.label = 'Stochastic';
workflowConfig.precompute.robustPlans.robust_3.objectiveSetName = 'Stochastic';
workflowConfig.precompute.robustPlans.robust_3.scenario.mode = 'wcScen';
workflowConfig.precompute.robustPlans.robust_3.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_3.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_3.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_3.scenario.wcSigma = 1.0;

workflowConfig.precompute.robustPlans.robust_4.label = 'c-Minimax';
workflowConfig.precompute.robustPlans.robust_4.objectiveSetName = 'cMinimax';
workflowConfig.precompute.robustPlans.robust_4.scenario.mode = 'wcScen';
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
workflowConfig.precompute.robustPlans.robust_5.scenario.mode = 'wcScen';
workflowConfig.precompute.robustPlans.robust_5.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_5.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_5.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_5.scenario.wcSigma = 1.0;

workflowConfig.precompute.robustPlans.robust_6.label = 'INTERVAL2';
workflowConfig.precompute.robustPlans.robust_6.objectiveSetName = 'Interval2';
workflowConfig.precompute.robustPlans.robust_6.scenario.mode = 'wcScen';
workflowConfig.precompute.robustPlans.robust_6.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_6.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_6.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_6.scenario.wcSigma = 1.0;
workflowConfig.precompute.robustPlans.robust_6.variants(1).id = 'theta_10';
workflowConfig.precompute.robustPlans.robust_6.variants(1).label = 'theta1=10';
workflowConfig.precompute.robustPlans.robust_6.variants(1).theta1 = 10;

workflowConfig.precompute.robustPlans.robust_7.label = 'INTERVAL3';
workflowConfig.precompute.robustPlans.robust_7.objectiveSetName = 'Interval3';
workflowConfig.precompute.robustPlans.robust_7.scenario.mode = 'wcScen';
workflowConfig.precompute.robustPlans.robust_7.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_7.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_7.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_7.scenario.wcSigma = 1.0;
workflowConfig.precompute.robustPlans.robust_7.variants(1).id = 'theta1_10_theta2_1';
workflowConfig.precompute.robustPlans.robust_7.variants(1).label = 'theta1=10 - theta2=1';
workflowConfig.precompute.robustPlans.robust_7.variants(1).theta1 = 10;
workflowConfig.precompute.robustPlans.robust_7.variants(1).theta2 = 1;

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
workflowConfig.sampling.sampling_scen_mode = 'wcScen';
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
