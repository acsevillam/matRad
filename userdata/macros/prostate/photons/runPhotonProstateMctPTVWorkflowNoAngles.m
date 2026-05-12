function runPhotonProstateMctPTVWorkflowNoAngles(varargin)
% runPhotonProstateMctPTVWorkflowNoAngles Run photon MCT prostate PTV workflow
% without gantry/couch angle perturbations in sampling.

macroFolder = fileparts(mfilename('fullpath'));
macroRoot = findMacroRoot(macroFolder);
helperFolder = fullfile(macroRoot,'helpers');
sharedFolder = fullfile(macroRoot,'shared','prostate');

if exist(helperFolder,'dir') == 7 && ~any(strcmp(strsplit(path,pathsep),helperFolder))
    addpath(helperFolder);
end
if exist(sharedFolder,'dir') ~= 7
    error('planWorkflow:macros:MissingSharedMacro', ...
        'Shared prostate macro folder not found: %s.',sharedFolder);
end
if ~any(strcmp(strsplit(path,pathsep),sharedFolder))
    addpath(sharedFolder);
end
if exist('ensurePlanWorkflowAvailable','file') == 2
    ensurePlanWorkflowAvailable('runPhotonProstateMctPTVWorkflowNoAngles');
end

userDataRoot = fileparts(macroRoot);
prepareConfig = photonProstatePrepareConfig( ...
    prostateMctPrepareConfig('PTV_001'));

macroDefaults = struct();
macroDefaults.caseID = fieldOrDefault(prepareConfig,'caseID','1_mct');
macroDefaults.rootPath = userDataRoot;
macroDefaults.cacheRootPath = fullfile(userDataRoot,'output','cache');
macroDefaults.randomSeed = 9999;
[macroArgs,randomSeed] = splitRandomSeedOption(macroDefaults.randomSeed,varargin{:});
macroOptions = planWorkflow.gui.PlanPresetWriter.parseMacroOptions( ...
    macroDefaults,macroArgs{:});

workflowConfig = struct();

workflowConfig.rootPath = macroOptions.rootPath;
workflowConfig.outputRootPath = fullfile(macroOptions.rootPath,'output');
workflowConfig.patientDataPath = fullfile(macroOptions.rootPath,'patients');
workflowConfig.cacheRootPath = macroOptions.cacheRootPath;

workflowConfig.prepare = prepareConfig;
workflowConfig.prepare.caseID = macroOptions.caseID;

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

workflowConfig.precompute.useCache = true;
workflowConfig.precompute.writeCache = true;

workflowConfig.pullDose.step1Enabled = true;
workflowConfig.pullDose.step1Target = {'CTV'};
workflowConfig.pullDose.step1Criteria = {'COV1'};
workflowConfig.pullDose.step1Limit = 0.9;
workflowConfig.pullDose.step1Start = 10;
workflowConfig.pullDose.step2Enabled = false;
workflowConfig.pullDose.maxIterations = 100;

workflowConfig.optimize.optimizer = 'IPOPT';

workflowConfig.sampling.sampling_linkToOptimization = true;
workflowConfig.sampling.sampling_scen_mode = 'random';
workflowConfig.sampling.sampling_ctActive = false;
workflowConfig.sampling.sampling_ctReferenceScenId = 1;
workflowConfig.sampling.sampling_setupActive = true;
workflowConfig.sampling.sampling_shiftSD = [5 5 5];
workflowConfig.sampling.sampling_wcSigma = 1.5;
workflowConfig.sampling.sampling_gantryActive = false;
workflowConfig.sampling.sampling_gantryAngleSD = 0;
workflowConfig.sampling.sampling_couchActive = false;
workflowConfig.sampling.sampling_couchAngleSD = 0;
workflowConfig.sampling.sampling_randomSeed = randomSeed;

workflowConfig.analysis.evaluationMode = 'total';
workflowConfig.analysis.gammaWindow = [0 1];
workflowConfig.analysis.gammaCriteria = [3 3];
workflowConfig.analysis.robustnessCriteria = [5 5];
workflowConfig.analysis.robustnessTargetMode = 'include';
workflowConfig.analysis.robustnessTargets = {'CTV'};

workflow = planWorkflow.Workflow(workflowConfig);
cleanupObj = onCleanup(@() workflow.releaseMemory()); %#ok<NASGU>

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

function value = fieldOrDefault(s,fieldName,defaultValue)

if isfield(s,fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end

end
