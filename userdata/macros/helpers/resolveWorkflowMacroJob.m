function jobPlan = resolveWorkflowMacroJob(job,varargin)
% resolveWorkflowMacroJob Resolve a serial macro job without running stages.

ensureJobDependencies('resolveWorkflowMacroJob');
job = normalizeJob(job);
runs = expandJobRuns(job);

jobPlan = struct();
jobPlan.id = job.id;
jobPlan.description = job.description;
jobPlan.stopOnError = job.stopOnError;
jobPlan.runs = repmat(emptyResolvedRun(),1,numel(runs));

for runIx = 1:numel(runs)
    runConfig = runs(runIx);
    runConfig = normalizeRunConfig(runConfig,runIx);
    specId = macroSpecId( ...
        runConfig.site,runConfig.particleType,runConfig.caseID, ...
        runConfig.robustness,runConfig.samplingProfile);
    macroSpecCatalog(specId,'profile',runConfig.profile);

    runOptions = mergeStructs(job.options,runConfig.options);
    args = [ ...
        {'profile',runConfig.profile, ...
         'openGui',runConfig.openGui, ...
         'optimizationScenario',runConfig.optimizationScenario, ...
         'samplingScenario',runConfig.samplingScenario}, ...
        structToNameValue(runOptions),varargin];

    jobPlan.runs(runIx).index = runIx;
    jobPlan.runs(runIx).label = runConfig.label;
    jobPlan.runs(runIx).specId = specId;
    jobPlan.runs(runIx).profile = runConfig.profile;
    jobPlan.runs(runIx).openGui = runConfig.openGui;
    jobPlan.runs(runIx).selector = struct( ...
        'site',runConfig.site, ...
        'particleType',runConfig.particleType, ...
        'caseID',runConfig.caseID, ...
        'robustness',runConfig.robustness, ...
        'samplingProfile',runConfig.samplingProfile);
    jobPlan.runs(runIx).optimizationScenario = ...
        runConfig.optimizationScenario;
    jobPlan.runs(runIx).samplingScenario = runConfig.samplingScenario;
    jobPlan.runs(runIx).options = runOptions;
    jobPlan.runs(runIx).args = args;
end

end

function ensureJobDependencies(context)

macroRoot = fileparts(fileparts(mfilename('fullpath')));
specFolder = fullfile(macroRoot,'shared','specs');
if exist(specFolder,'dir') ~= 7
    error('planWorkflow:macros:MissingSpecCatalog', ...
        'Macro spec catalog folder not found: %s.',specFolder);
end
if ~any(strcmp(strsplit(path,pathsep),specFolder))
    addpath(specFolder);
end
ensurePlanWorkflowAvailable(context);

end

function job = normalizeJob(job)

if ~isstruct(job) || ~isscalar(job)
    error('planWorkflow:macros:InvalidMacroJob', ...
        'Macro job must be a scalar struct.');
end
assertKnownFields(job,allowedJobFields(),'job');

job = normalizeJobHeader(job);
job = normalizeJobExecutionDefaults(job);
job = normalizeJobScenarioDefaults(job);
job = normalizeJobOptions(job);

[hasRuns,hasBases] = macroJobShape(job);
validateMacroJobShape(hasRuns,hasBases);
job = normalizeMacroJobDefinitions(job,hasRuns);

end

function job = normalizeJobHeader(job)

if ~isfield(job,'id') || isempty(job.id)
    job.id = 'macroJob';
end
job.id = normalizeText(job.id,'job.id');

if ~isfield(job,'description') || isempty(job.description)
    job.description = '';
end
job.description = normalizeOptionalText( ...
    job.description,'job.description');

end

function job = normalizeJobExecutionDefaults(job)

if ~isfield(job,'profile') || isempty(job.profile)
    job.profile = 'prod';
end
job.profile = normalizeProfile(job.profile);

if ~isfield(job,'openGui') || isempty(job.openGui)
    job.openGui = false;
end
job.openGui = planWorkflow.config.ConfigValue.logicalScalar( ...
    job.openGui,'job.openGui', ...
    'planWorkflow:macros:InvalidMacroJob');

if ~isfield(job,'stopOnError') || isempty(job.stopOnError)
    job.stopOnError = true;
end
job.stopOnError = planWorkflow.config.ConfigValue.logicalScalar( ...
    job.stopOnError,'job.stopOnError', ...
    'planWorkflow:macros:InvalidMacroJob');

end

function job = normalizeJobScenarioDefaults(job)

if ~isfield(job,'samplingProfile') || isempty(job.samplingProfile)
    job.samplingProfile = 'default';
end
job.samplingProfile = normalizeText( ...
    job.samplingProfile,'job.samplingProfile');

if ~isfield(job,'optimizationScenario') || ...
        isempty(job.optimizationScenario)
    job.optimizationScenario = struct();
end
job.optimizationScenario = normalizeScalarStruct( ...
    job.optimizationScenario,'job.optimizationScenario');

if ~isfield(job,'samplingScenario') || isempty(job.samplingScenario)
    job.samplingScenario = struct();
end
job.samplingScenario = normalizeScalarStruct( ...
    job.samplingScenario,'job.samplingScenario');

end

function job = normalizeJobOptions(job)

if ~isfield(job,'options') || isempty(job.options)
    job.options = struct();
end
job.options = normalizeScalarStruct(job.options,'job.options');

end

function [hasRuns,hasBases] = macroJobShape(job)

hasRuns = isfield(job,'runs') && ~isempty(job.runs);
hasBases = hasTopLevelMacroBase(job) || ...
    hasNonemptyField(job,'base') || hasNonemptyField(job,'bases');

end

function tf = hasTopLevelMacroBase(job)

tf = isfield(job,'site') || isfield(job,'particleType') || ...
    isfield(job,'caseID') || isfield(job,'robustness');

end

function tf = hasNonemptyField(value,fieldName)

tf = isfield(value,fieldName) && ~isempty(value.(fieldName));

end

function validateMacroJobShape(hasRuns,hasBases)

if hasRuns && hasBases
    error('planWorkflow:macros:InvalidMacroJob', ...
        'Macro job must define either runs or base/bases, not both.');
end
if ~hasRuns && ~hasBases
    error('planWorkflow:macros:InvalidMacroJob', ...
        ['Macro job must define either runs or base/bases with ' ...
         'parameterSets.']);
end

end

function job = normalizeMacroJobDefinitions(job,hasRuns)

if hasRuns
    if ~isstruct(job.runs)
        error('planWorkflow:macros:InvalidMacroJob', ...
            'job.runs must be a non-empty struct array.');
    end
    job.runs = normalizePatchArray(job.runs,'job.runs');
else
    if isfield(job,'base') && ~isempty(job.base)
        job.base = normalizePatchArray(job.base,'job.base');
    end
    if isfield(job,'bases') && ~isempty(job.bases)
        job.bases = normalizePatchArray(job.bases,'job.bases');
    end
    if ~isfield(job,'parameterSets') || isempty(job.parameterSets)
        job.parameterSets = struct();
    else
        job.parameterSets = normalizePatchArray( ...
            job.parameterSets,'job.parameterSets');
    end
end

end

function runs = expandJobRuns(job)

defaults = jobDefaults(job);
if isfield(job,'runs') && ~isempty(job.runs)
    runs = repmat(defaults,1,numel(job.runs));
    for runIx = 1:numel(job.runs)
        runs(runIx) = mergeStructs(defaults,job.runs(runIx));
    end
    return;
end

bases = collectBases(job);
parameterSets = collectParameterSets(job);
runs = repmat(defaults,1,numel(bases)*numel(parameterSets));
outIx = 1;
for baseIx = 1:numel(bases)
    baseConfig = mergeStructs(defaults,bases{baseIx});
    for parameterSetIx = 1:numel(parameterSets)
        runConfig = mergeStructs( ...
            baseConfig,parameterSets(parameterSetIx));
        if ~isfield(parameterSets(parameterSetIx),'label') || ...
                isempty(parameterSets(parameterSetIx).label)
            runConfig.label = defaultRunLabel( ...
                baseConfig,parameterSets(parameterSetIx), ...
                baseIx,parameterSetIx);
        end
        runs(outIx) = runConfig;
        outIx = outIx + 1;
    end
end

end

function bases = collectBases(job)

baseList = {};
if isfield(job,'base') && ~isempty(job.base)
    for baseIx = 1:numel(job.base)
        baseList{end + 1} = job.base(baseIx); %#ok<AGROW>
    end
end
if isfield(job,'bases') && ~isempty(job.bases)
    for baseIx = 1:numel(job.bases)
        baseList{end + 1} = job.bases(baseIx); %#ok<AGROW>
    end
end
if isempty(baseList)
    baseList = {struct()};
end

bases = baseList;

end

function parameterSets = collectParameterSets(job)

if isfield(job,'parameterSets') && ~isempty(job.parameterSets)
    parameterSets = job.parameterSets;
else
    parameterSets = struct();
end

end

function label = defaultRunLabel(baseConfig,parameterSet,baseIx, ...
        parameterSetIx)

if isfield(parameterSet,'label') && ~isempty(parameterSet.label)
    label = parameterSet.label;
elseif isfield(parameterSet,'robustness') && ...
        ~isempty(parameterSet.robustness)
    label = parameterSet.robustness;
elseif isfield(parameterSet,'caseID') && ~isempty(parameterSet.caseID)
    label = parameterSet.caseID;
elseif isfield(baseConfig,'robustness') && ~isempty(baseConfig.robustness)
    label = baseConfig.robustness;
else
    label = sprintf('base%d_parameterSet%d',baseIx,parameterSetIx);
end
label = char(string(label));
if isfield(baseConfig,'label') && ~isempty(baseConfig.label)
    baseLabel = char(string(baseConfig.label));
    if ~strcmp(baseLabel,label)
        label = [baseLabel '_' label];
    end
end

end

function patches = normalizePatchArray(patches,context)

if ~isstruct(patches)
    error('planWorkflow:macros:InvalidMacroJob', ...
        '%s must be a struct array.',context);
end
if ~isfield(patches,'options')
    for patchIx = 1:numel(patches)
        patches(patchIx).options = struct();
    end
end
for patchIx = 1:numel(patches)
    patches(patchIx) = normalizeRunPatch( ...
        patches(patchIx),sprintf('%s(%d)',context,patchIx));
end

end

function patch = normalizeRunPatch(patch,context)

assertKnownFields(patch,allowedRunFields(),context);
if ~isfield(patch,'options') || isempty(patch.options)
    patch.options = struct();
end
patch.options = normalizeScalarStruct( ...
    patch.options,sprintf('%s.options',context));

end

function defaults = jobDefaults(job)

defaults = struct();
defaults.label = '';
defaults.profile = job.profile;
defaults.openGui = job.openGui;
defaults.site = optionalDefault(job,'site','');
defaults.particleType = optionalDefault(job,'particleType','');
defaults.caseID = optionalDefault(job,'caseID','');
defaults.robustness = optionalDefault(job,'robustness','');
defaults.samplingProfile = job.samplingProfile;
defaults.optimizationScenario = job.optimizationScenario;
defaults.samplingScenario = job.samplingScenario;
defaults.options = struct();

end

function value = optionalDefault(job,fieldName,defaultValue)

if isfield(job,fieldName) && ~isempty(job.(fieldName))
    value = job.(fieldName);
else
    value = defaultValue;
end

end

function runConfig = normalizeRunConfig(runConfig,runIx)

runConfig.profile = normalizeProfile(runConfig.profile);
runConfig.openGui = planWorkflow.config.ConfigValue.logicalScalar( ...
    runConfig.openGui,sprintf('job.runs(%d).openGui',runIx), ...
    'planWorkflow:macros:InvalidMacroJob');
runConfig.site = normalizeText( ...
    runConfig.site,sprintf('job.runs(%d).site',runIx));
runConfig.particleType = normalizeText( ...
    runConfig.particleType,sprintf('job.runs(%d).particleType',runIx));
runConfig.caseID = normalizeText( ...
    runConfig.caseID,sprintf('job.runs(%d).caseID',runIx));
runConfig.robustness = normalizeText( ...
    runConfig.robustness,sprintf('job.runs(%d).robustness',runIx));
runConfig.samplingProfile = normalizeText( ...
    runConfig.samplingProfile, ...
    sprintf('job.runs(%d).samplingProfile',runIx));
runConfig.optimizationScenario = normalizeScalarStruct( ...
    runConfig.optimizationScenario, ...
    sprintf('job.runs(%d).optimizationScenario',runIx));
runConfig.samplingScenario = normalizeScalarStruct( ...
    runConfig.samplingScenario, ...
    sprintf('job.runs(%d).samplingScenario',runIx));
runConfig.options = normalizeScalarStruct( ...
    runConfig.options,sprintf('job.runs(%d).options',runIx));

if isempty(runConfig.label)
    runConfig.label = runConfig.robustness;
end
runConfig.label = normalizeText( ...
    runConfig.label,sprintf('job.runs(%d).label',runIx));

end

function fields = allowedJobFields()

fields = {'id','description','profile','openGui','stopOnError', ...
    'site','particleType','caseID','robustness','samplingProfile', ...
    'optimizationScenario','samplingScenario','options','base','bases', ...
    'parameterSets','runs'};

end

function fields = allowedRunFields()

fields = {'label','profile','openGui','site','particleType','caseID', ...
    'robustness','samplingProfile','optimizationScenario', ...
    'samplingScenario','options'};

end

function assertKnownFields(value,allowed,context)

fields = fieldnames(value);
for fieldIx = 1:numel(fields)
    if ~any(strcmp(fields{fieldIx},allowed))
        error('planWorkflow:macros:InvalidMacroJob', ...
            'Unknown field "%s" in %s.',fields{fieldIx},context);
    end
end

end

function value = normalizeText(value,context)

if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('planWorkflow:macros:InvalidMacroJob', ...
        '%s must be a non-empty text scalar.',context);
end
value = char(string(value));
if isempty(strtrim(value))
    error('planWorkflow:macros:InvalidMacroJob', ...
        '%s must be a non-empty text scalar.',context);
end

end

function value = normalizeOptionalText(value,context)

if isempty(value)
    value = '';
    return;
end
value = normalizeText(value,context);

end

function profile = normalizeProfile(profile)

profile = lower(normalizeText(profile,'profile'));
if ~any(strcmp(profile,{'prod','testing'}))
    error('planWorkflow:macros:InvalidMacroProfile', ...
        'Macro profile must be "prod" or "testing".');
end

end

function value = normalizeScalarStruct(value,context)

if isempty(value)
    value = struct();
end
if ~isstruct(value) || ~isscalar(value)
    error('planWorkflow:macros:InvalidMacroJob', ...
        '%s must be a scalar struct.',context);
end

end

function values = structToNameValue(options)

fields = fieldnames(options);
values = cell(1,2*numel(fields));
outIx = 1;
for fieldIx = 1:numel(fields)
    values{outIx} = fields{fieldIx};
    values{outIx + 1} = options.(fields{fieldIx});
    outIx = outIx + 2;
end

end

function merged = mergeStructs(base,patch)

if nargin < 1 || isempty(base)
    base = struct();
end
if nargin < 2 || isempty(patch)
    merged = base;
    return;
end
if ~isstruct(base) || ~isstruct(patch) || ...
        ~(isscalar(base) && isscalar(patch))
    merged = patch;
    return;
end

merged = base;
fields = fieldnames(patch);
for fieldIx = 1:numel(fields)
    fieldName = fields{fieldIx};
    if isfield(merged,fieldName)
        merged.(fieldName) = mergeStructs( ...
            merged.(fieldName),patch.(fieldName));
    else
        merged.(fieldName) = patch.(fieldName);
    end
end

end

function run = emptyResolvedRun()

run = struct();
run.index = [];
run.label = '';
run.specId = '';
run.profile = '';
run.openGui = false;
run.selector = struct();
run.optimizationScenario = struct();
run.samplingScenario = struct();
run.options = struct();
run.args = {};

end
