function runProstateCOWCWorkflowCore(prepareConfig,varargin)
% runProstateCOWCWorkflowCore Run the shared prostate COWC workflow.
% Usage: runProstateCOWCWorkflowCore(prepareConfig,'caseID','1_mct','rootPath',userDataRoot,'cacheRootPath',fullfile(userDataRoot,'output','cache'),'randomSeed',[])

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
    ensurePlanWorkflowAvailable('runProstateCOWCWorkflowCore');
end
userDataRoot = fileparts(macroRoot);

macroDefaults = struct();
macroDefaults.caseID = fieldOrDefault(prepareConfig,'caseID','1_mct');
macroDefaults.rootPath = userDataRoot;
macroDefaults.cacheRootPath = fullfile(userDataRoot,'output','cache');
macroDefaults.randomSeed = 9999;
[macroArgs,randomSeed] = splitRandomSeedOption(macroDefaults.randomSeed,varargin{:});
macroOptions = planWorkflow.gui.PlanPresetWriter.parseMacroOptions( ...
    macroDefaults,macroArgs{:});

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

workflowConfig.precompute.robustPlans.robust_1.label = 'Minimax';
workflowConfig.precompute.robustPlans.robust_1.objectiveSetName = 'Minimax';
workflowConfig.precompute.robustPlans.robust_1.scenario.mode = 'wcScen';
workflowConfig.precompute.robustPlans.robust_1.scenario.ctActive = false;
workflowConfig.precompute.robustPlans.robust_1.scenario.ctReferenceScenId = 1;
workflowConfig.precompute.robustPlans.robust_1.scenario.setupActive = true;
workflowConfig.precompute.robustPlans.robust_1.scenario.shiftSD = [5 5 5];
workflowConfig.precompute.robustPlans.robust_1.scenario.rangeActive = true;
workflowConfig.precompute.robustPlans.robust_1.scenario.rangeAbsSD = 1;
workflowConfig.precompute.robustPlans.robust_1.scenario.rangeRelSD = 3.5;
workflowConfig.precompute.robustPlans.robust_1.scenario.numOfRangeGridPoints = 3;

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
workflowConfig.sampling.sampling_scen_mode = 'random';
workflowConfig.sampling.sampling_ctActive = false;
workflowConfig.sampling.sampling_ctReferenceScenId = 1;
workflowConfig.sampling.sampling_setupActive = true;
workflowConfig.sampling.sampling_shiftSD = [5 5 5];
workflowConfig.sampling.sampling_wcSigma = 1.5;
workflowConfig.sampling.sampling_rangeActive = true;
workflowConfig.sampling.sampling_rangeAbsSD = 1;
workflowConfig.sampling.sampling_rangeRelSD = 3.5;
workflowConfig.sampling.sampling_gantryActive = true;
workflowConfig.sampling.sampling_gantryAngleSD = 1;
workflowConfig.sampling.sampling_couchActive = true;
workflowConfig.sampling.sampling_couchAngleSD = 1;
workflowConfig.sampling.sampling_randomSeed = randomSeed;

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

function [macroArgs,randomSeed] = splitRandomSeedOption(defaultValue,varargin)

macroArgs = varargin;
randomSeed = defaultValue;
if isempty(varargin)
    return;
end

if numel(varargin) == 1 && isstruct(varargin{1})
    options = varargin{1};
    if isfield(options,'randomSeed')
        randomSeed = options.randomSeed;
        options = rmfield(options,'randomSeed');
    end
    macroArgs = {options};
    return;
end

if ~isNameValueArgs(varargin)
    return;
end

macroArgs = {};
i = 1;
while i <= numel(varargin)
    name = char(string(varargin{i}));
    value = varargin{i + 1};
    if strcmp(name,'randomSeed')
        randomSeed = value;
    else
        macroArgs = [macroArgs varargin(i:i+1)]; %#ok<AGROW>
    end
    i = i + 2;
end

end

function tf = isNameValueArgs(args)

tf = false;
if mod(numel(args),2) ~= 0
    return;
end
if isempty(args)
    return;
end
firstName = args{1};
if ~(ischar(firstName) || (isstring(firstName) && isscalar(firstName)))
    return;
end
knownNames = {'caseID','rootPath','cacheRootPath','randomSeed'};
tf = any(strcmp(char(string(firstName)),knownNames));

end
