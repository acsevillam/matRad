function result = runBreastPhotonMctMultipleCtWorkflow(varargin)
% runBreastPhotonMctMultipleCtWorkflow Configure and run a CT-aware breast MacroSpec.

% Path setup.
macroRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(macroRoot,'helpers'));
specFolder = fullfile(macroRoot,'shared','specs');
if ~any(strcmp(strsplit(path,pathsep),specFolder))
    addpath(specFolder);
end
ensurePlanWorkflowAvailable('runBreastPhotonMctMultipleCtWorkflow');

% Execution defaults.
% profile: prod uses production defaults; testing applies explicit overrides.
% openGui: true opens the planWorkflow GUI before running stages.
profile = 'prod';
openGui = true;

% MacroSpec selectors.
% site: anatomy/catalog key (breast, prostate, head_and_neck).
% particleType: radiation mode used by prepare.radiationMode.
% caseID: patient/case identifier resolved under userdata.
% robustness: canonical robust plan key or multiple-plan selector.
% samplingProfile: sampling profile (default, noAngles).
site = 'breast';
particleType = 'photons';
caseID = '4136_mct';
robustness = 'multiple';
samplingProfile = 'default';

% Optimization scenario toggles.
% Controls robust scenarios used for precompute and optimization plans.
optimizationScenario = struct( ...
    'ctActive', false, ...
    'setupActive', true, ...
    'rangeActive', false, ...
    'gantryActive', false, ...
    'couchActive', false);
ctOptimizationScenario = optimizationScenario;
ctOptimizationScenario.ctActive = true;

% Sampling scenario toggles.
% Controls active uncertainty dimensions used by the sampling stage.
samplingScenario = struct( ...
    'ctActive', true, ...
    'setupActive', true, ...
    'rangeActive', false, ...
    'gantryActive', false, ...
    'couchActive', false);

% Resolve the MacroSpec and run the complete workflow.
specId = macroSpecId( ...
    site,particleType,caseID,robustness,samplingProfile);
[profile,args] = parseMacroProfileOption(profile,varargin{:});
robustPlans = robustPlansWithPlanSpecificCt( ...
    specId,profile,optimizationScenario,ctOptimizationScenario);
result = runWorkflowMacroSpec( ...
    specId,'profile',profile,'openGui',openGui, ...
    'optimizationScenario',optimizationScenario, ...
    'samplingScenario',samplingScenario, ...
    'allowCustomRobustPlans',true, ...
    'precompute',struct('robustPlans',robustPlans),args{:});
end

function robustPlans = robustPlansWithPlanSpecificCt( ...
        specId,profile,optimizationScenario,ctOptimizationScenario)

spec = macroSpecCatalog(specId,'profile',profile);
robustPlans = planWorkflow.config.RobustPlanCatalog.select( ...
    spec.description,spec.planTemplate,spec.planKeys, ...
    'nominalScenario',spec.nominalScenario, ...
    'robustScenario',optimizationScenario, ...
    'radiationMode',spec.modality);

targetObjectiveSets = {'MeanVariance','Interval2','Interval3'};
matchedObjectiveSets = cell(1,0);
for planIx = 1:numel(robustPlans)
    objectiveSetName = char(robustPlans(planIx).objectiveSetName);
    if any(strcmp(objectiveSetName,targetObjectiveSets))
        robustPlans(planIx).scenario = ...
            planWorkflow.config.ScenarioSpec.normalize( ...
            ctOptimizationScenario,robustPlans(planIx).scenario, ...
            sprintf('optimizationScenario.%s',objectiveSetName));
        matchedObjectiveSets{end + 1} = objectiveSetName; %#ok<AGROW>
    end
end

missingObjectiveSets = setdiff( ...
    targetObjectiveSets,matchedObjectiveSets,'stable');
if ~isempty(missingObjectiveSets)
    error('planWorkflow:macros:MissingCtRobustPlans', ...
        'Missing expected CT optimization plan(s): %s.', ...
        strjoin(missingObjectiveSets,', '));
end
end
