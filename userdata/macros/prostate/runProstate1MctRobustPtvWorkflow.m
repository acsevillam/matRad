function runProstate1MctRobustPtvWorkflow(varargin)
% runProstate1MctRobustPtvWorkflow Run the 1_mct prostate robust PTV workflow.
% Usage: runProstate1MctRobustPtvWorkflow('rootPath',userDataRoot,'cacheRootPath',fullfile(userDataRoot,'output','cache'))

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
    ensurePlanWorkflowAvailable('runProstate1MctRobustPtvWorkflow');
end
userDataRoot = fileparts(macroRoot);

macroDefaults = struct();
macroDefaults.caseID = '1_mct';
macroDefaults.rootPath = userDataRoot;
macroDefaults.cacheRootPath = fullfile(userDataRoot,'output','cache');
macroDefaults.randomSeed = 9999;
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
workflowConfig.prepare.plan_template = 'PTV_001';
workflowConfig.prepare.radiationMode = 'photons';
workflowConfig.prepare.machine = 'Generic';
workflowConfig.prepare.bioModel = 'none';
workflowConfig.prepare.plan_beams = '9F';
workflowConfig.prepare.dicomMetadata.ctSeriesUIDs = { ...
    '1.2.246.352.221.4862331015009694712.7895934050572435376', ...
    '1.2.246.352.221.5242610798239234972.16653089788686535091', ...
    '1.2.246.352.221.4884747965284749950.2674045711264065698', ...
    '1.2.246.352.221.4763733919031700794.11500848528649628325', ...
    '1.2.246.352.221.5463932454974956859.8622902214826974644', ...
    '1.2.246.352.221.5512918193748952617.16453683212378280614'};
workflowConfig.prepare.dicomMetadata.rtssUIDs = { ...
    '1.2.246.352.221.5753679130948318511.11670799102824828582', ...
    '1.2.246.352.221.5709118728753457411.6048825263865128877', ...
    '1.2.246.352.221.5073606025232248792.2859148916510586007', ...
    '1.2.246.352.221.4871210827866228606.2203101012107793342', ...
    '1.2.246.352.221.4946890673530045409.2568569967401124497', ...
    '1.2.246.352.221.4705110154949056430.3657753086562001833'};
workflowConfig.prepare.resolution = [3 3 3];
workflowConfig.prepare.n_cores = 8;

% Precompute stage.
workflowConfig.precompute.doseResolution = [3 3 3];
workflowConfig.precompute.reference.label = 'Nominal';
workflowConfig.precompute.reference.strategy = 'none';
workflowConfig.precompute.reference.scenario.mode = 'nomScen';
workflowConfig.precompute.reference.scenario.ctActive = false;
workflowConfig.precompute.reference.scenario.ctReferenceScenId = 1;

workflowConfig.precompute.robustPlans.robust_1.label = 'Robust PTV';
workflowConfig.precompute.robustPlans.robust_1.objectiveSetName = 'robust_1';
workflowConfig.precompute.robustPlans.robust_1.strategy = 'COWC';
workflowConfig.precompute.robustPlans.robust_1.scenario.mode = 'wcScen';
workflowConfig.precompute.robustPlans.robust_1.scenario.ctActive = true;
workflowConfig.precompute.robustPlans.robust_1.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_1.scenario.shiftSD = [5 10 5];
workflowConfig.precompute.robustPlans.robust_1.scenario.wcSigma = 1.0;

workflowConfig.precompute.useCache = true;
workflowConfig.precompute.writeCache = true;

% Dose-pulling stage.
workflowConfig.pullDose.step1Enabled = false;
workflowConfig.pullDose.step2Enabled = false;
workflowConfig.pullDose.maxIterations = 100;
workflowConfig.pullDose.scale_factor = 1;

% Optimize stage.
workflowConfig.optimize.optimizer = 'IPOPT';

% Sample stage.
workflowConfig.sampling.linkToOptimization = true;
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
workflowConfig.analysis.robustnessTargets = {'CTV','PTV'};

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
