function spec = macroSpecCatalog(specId,varargin)
% macroSpecCatalog Resolve local workflow macro specs by ID.

if nargin < 1 || isempty(specId)
    spec = allSpecIds();
    return;
end
specId = char(string(specId));
if strcmp(specId,'all') || strcmp(specId,'ids')
    spec = allSpecIds();
    return;
end
if ~any(strcmp(specId,allSpecIds()))
    error('planWorkflow:macros:UnknownMacroSpec', ...
        'Unknown macro spec "%s".',specId);
end
[profile,unusedArgs] = parseMacroProfileOption('prod',varargin{:});
if ~isempty(unusedArgs)
    error('planWorkflow:macros:InvalidMacroSpecOptions', ...
        'macroSpecCatalog only accepts the "profile" option.');
end

parts = splitSpecId(specId);
site = parts{1};
modality = parts{2};
caseID = parts{3};
planAlias = parts{4};

spec = baseSpec(profile,site,modality,caseID,planAlias);
spec = planWorkflow.macros.MacroSpec.normalize(spec);

end

function ids = allSpecIds()

ids = {};
ids = [ids specIds('breast','photons','4136_mct', ...
    {'PTV','COWC','STOCH','cCOWC','PROB2','INTERVAL2', ...
    'INTERVAL3','multiple'})];
ids = [ids specIds('breast','photons','4136', ...
    {'multiple'})];
ids = [ids specIds('prostate','photons','1_mct', ...
    {'PTV','PTV_noAngles'})];
ids = [ids specIds('prostate','photons','3482', ...
    {'INTERVAL2','multiple'})];
ids = [ids specIds('prostate','protons','1_mct', ...
    {'PTV','COWC','COWC_noAngles'})];
ids = [ids specIds('prostate','protons','3482', ...
    {'multiple'})];
ids = [ids specIds('prostate','carbon','1_mct', ...
    {'PTV','COWC','COWC_noAngles'})];
ids = [ids specIds('prostate','helium','1_mct', ...
    {'PTV','COWC','COWC_noAngles'})];
ids = [ids specIds('head_and_neck','photons','2', ...
    {'PTV','INTERVAL2','INTERVAL3','multiple'})];

end

function ids = specIds(site,modality,caseID,planAliases)

ids = cell(1,numel(planAliases));
for aliasIx = 1:numel(planAliases)
    ids{aliasIx} = strjoin( ...
        {site,modality,caseID,planAliases{aliasIx}},'.');
end

end

function parts = splitSpecId(specId)

parts = strsplit(specId,'.');
if numel(parts) ~= 4
    error('planWorkflow:macros:InvalidMacroSpecId', ...
        ['Macro spec IDs must use ' ...
         '<site>.<modality>.<caseID>.<plan>.']);
end

end

function spec = applyProfileOverrides(spec,profile)

switch profile
    case 'prod'
        return;
    case 'testing'
        spec.precompute.doseResolution = [5 5 5];
    otherwise
        error('planWorkflow:macros:InvalidMacroProfile', ...
            'Macro profile must be "prod" or "testing".');
end

end

function spec = baseSpec(profile,site,modality,caseID,planAlias)

[profile,unusedArgs] = parseMacroProfileOption(profile);
if ~isempty(unusedArgs)
    error('planWorkflow:macros:InvalidMacroProfile', ...
        'Macro profile must be "prod" or "testing".');
end
planConfig = planAliasConfig(site,planAlias);
spec = struct();
spec.baseId = strjoin({site,modality,caseID,planAlias},'.');
spec.id = strjoin({profile,spec.baseId},'.');
spec.profile = profile;
spec.site = site;
spec.description = descriptionForSite(site);
spec.modality = modality;
spec.caseID = caseID;
spec.planTemplate = planConfig.template;
spec.planKeys = planConfig.planKeys;
spec.prepare = prepareConfig(site,modality,caseID,planConfig.template);
spec.nominalScenario = nominalScenario();
spec.robustScenario = robustScenario(site,planAlias);
spec.precompute = precomputeDefaults();
spec.reference = referenceDefaults();
spec.pullDose = pullDoseDefaults();
spec.optimize = struct('optimizer','IPOPT');
spec.sampling = samplingDefaults(site,planAlias);
spec.analysis = analysisDefaults();
spec.executionMode = 'run';
spec.openGui = true;
spec.lockPlanSet = true;
spec.allowCustomRobustPlans = false;
spec.randomSeed = randomSeedForPlan(site,planAlias);
spec.rootPath = userDataRoot();
spec = applyProfileOverrides(spec,profile);

end

function description = descriptionForSite(site)

switch site
    case 'breast'
        description = 'breast';
    case 'prostate'
        description = 'prostate';
    case 'head_and_neck'
        description = 'h&n';
    otherwise
        error('planWorkflow:macros:InvalidMacroSite', ...
            'Unsupported macro site "%s".',site);
end

end

function config = planAliasConfig(site,planAlias)

switch planAlias
    case {'PTV','PTV_noAngles'}
        config.template = 'PTV_001';
        config.planKeys = {'PTV'};
    case {'COWC','COWC_noAngles'}
        config.template = 'COWC_001';
        config.planKeys = {'COWC'};
    case 'STOCH'
        config.template = 'STOCH_001';
        config.planKeys = {'STOCH'};
    case 'cCOWC'
        config.template = 'cCOWC_001';
        config.planKeys = {'cCOWC'};
    case 'PROB2'
        config.template = 'PROB2_001';
        config.planKeys = {'PROB2'};
    case 'INTERVAL2'
        config.template = 'interval2_001';
        config.planKeys = {'INTERVAL2'};
    case 'INTERVAL3'
        config.template = 'interval3_001';
        config.planKeys = {'INTERVAL3'};
    case 'multiple'
        config.template = 'comparison_001';
        config.planKeys = {'all'};
    otherwise
        error('planWorkflow:macros:InvalidMacroPlan', ...
            'Unsupported macro plan "%s".',planAlias);
end

if strcmp(site,'head_and_neck') && ...
        ~any(strcmp(planAlias,{'PTV','INTERVAL2','INTERVAL3','multiple'}))
    error('planWorkflow:macros:InvalidMacroPlan', ...
        'Unsupported H&N macro plan "%s".',planAlias);
end

end

function prepare = prepareConfig(site,modality,caseID,planTemplate)

prepare = struct();
prepare.caseID = caseID;
prepare.AcquisitionType = 'dicom';
prepare.hlutFileName = 'matRad_default.hlut';
prepare.description = descriptionForSite(site);
prepare.plan_template = planTemplate;
prepare.radiationMode = modality;
prepare.machine = machineForModality(modality);
prepare.bioModel = bioModelForModality(modality);
prepare.quantityOpt = quantityForModality(modality);
prepare.plan_beams = beamSetForSiteModality(site,modality);
prepare.dicomMetadata = dicomMetadata(site,caseID);
prepare.resolution = [3 3 3];

end

function machine = machineForModality(modality)

switch modality
    case {'photons','protons','carbon','helium'}
        machine = 'Generic';
    otherwise
        error('planWorkflow:macros:InvalidMacroModality', ...
            'Unsupported macro modality "%s".',modality);
end

end

function bioModel = bioModelForModality(modality)

switch modality
    case 'photons'
        bioModel = 'none';
    case 'protons'
        bioModel = 'constRBE';
    case 'carbon'
        bioModel = 'LEM';
    case 'helium'
        bioModel = 'HEL';
end

end

function quantity = quantityForModality(modality)

if strcmp(modality,'photons')
    quantity = 'physicalDose';
else
    quantity = 'RBExD';
end

end

function beamSet = beamSetForSiteModality(site,modality)

if any(strcmp(modality,{'protons','carbon','helium'}))
    beamSet = '2F';
    return;
end
switch site
    case 'breast'
        beamSet = '7F';
    case 'prostate'
        beamSet = '9F';
    case 'head_and_neck'
        beamSet = '7F';
end

end

function metadata = dicomMetadata(site,caseID)

metadata = struct();
if strcmp(site,'breast') && strcmp(caseID,'4136_mct')
    metadata.ctSeriesUIDs = { ...
        '1.2.752.243.1.1.20250926114032925.4100.14786', ...
        '1.2.752.243.1.1.20250926114142288.1030.41286', ...
        '1.2.752.243.1.1.20250926114610558.2260.80247', ...
        '1.2.752.243.1.1.20250926114701357.2890.12431', ...
        '1.2.752.243.1.1.20250926114759405.3520.62654'};
    metadata.rtssUID = ...
        '1.2.752.243.1.1.20250926114032968.9900.45533';
elseif strcmp(site,'prostate') && strcmp(caseID,'1_mct')
    metadata.ctSeriesUIDs = { ...
        '1.2.246.352.221.4862331015009694712.7895934050572435376', ...
        '1.2.246.352.221.5242610798239234972.16653089788686535091', ...
        '1.2.246.352.221.4884747965284749950.2674045711264065698', ...
        '1.2.246.352.221.4763733919031700794.11500848528649628325', ...
        '1.2.246.352.221.5463932454974956859.8622902214826974644', ...
        '1.2.246.352.221.5512918193748952617.16453683212378280614'};
    metadata.rtssUIDs = { ...
        '1.2.246.352.221.5753679130948318511.11670799102824828582', ...
        '1.2.246.352.221.5709118728753457411.6048825263865128877', ...
        '1.2.246.352.221.5073606025232248792.2859148916510586007', ...
        '1.2.246.352.221.4871210827866228606.2203101012107793342', ...
        '1.2.246.352.221.4946890673530045409.2568569967401124497', ...
        '1.2.246.352.221.4705110154949056430.3657753086562001833'};
end

end

function scenario = nominalScenario()

scenario = planWorkflow.config.ScenarioSpec.defaults('nomScen');
scenario.ctActive = false;
scenario.ctReferenceScenId = 1;
scenario.setupActive = false;
scenario.rangeActive = false;
scenario.gantryActive = false;
scenario.couchActive = false;

end

function scenario = robustScenario(site,planAlias)

if strcmp(site,'prostate') && contains(planAlias,'COWC')
    scenario = planWorkflow.config.ScenarioSpec.defaults('wcScen');
    scenario.ctActive = false;
    scenario.ctReferenceScenId = 1;
    scenario.setupActive = true;
    scenario.shiftSD = [5 5 5];
    scenario.rangeActive = true;
    scenario.rangeAbsSD = 1;
    scenario.rangeRelSD = 3.5;
    scenario.numOfRangeGridPoints = 3;
    scenario.gantryActive = false;
    scenario.couchActive = false;
    return;
end

scenario = planWorkflow.config.ScenarioSpec.defaults('impScen5');
scenario.ctActive = strcmp(site,'prostate') || strcmp(site,'head_and_neck');
scenario.ctReferenceScenId = 1;
scenario.setupActive = true;
scenario.rangeActive = false;
scenario.gantryActive = false;
scenario.couchActive = false;
switch site
    case 'breast'
        scenario.shiftSD = [4 8 6];
    case 'prostate'
        scenario.shiftSD = [5 5 5];
    case 'head_and_neck'
        scenario.shiftSD = [3 3 3];
end
scenario.wcSigma = 1.0;

end

function precompute = precomputeDefaults()

precompute = struct();
precompute.doseResolution = [3 3 3];
precompute.useCache = true;
precompute.writeCache = true;

end

function reference = referenceDefaults()

reference = struct();
reference.label = 'Nominal';
reference.scenario = nominalScenario();

end

function pullDose = pullDoseDefaults()

pullDose = struct();
pullDose.step1Enabled = true;
pullDose.step1Target = {'CTV'};
pullDose.step1Criteria = {'COV1'};
pullDose.step1Limit = 0.9;
pullDose.step1Start = 10;
pullDose.step2Enabled = false;
pullDose.maxIterations = 100;

end

function sampling = samplingDefaults(site,planAlias)

sampling = struct();
sampling.sampling_linkToOptimization = true;

if strcmp(site,'prostate') && strcmp(planAlias,'PTV')
    sampling = prostateAngularSampling(true,false);
elseif strcmp(site,'prostate') && strcmp(planAlias,'PTV_noAngles')
    sampling = prostateAngularSampling(false,false);
elseif strcmp(site,'prostate') && strcmp(planAlias,'COWC')
    sampling = prostateAngularSampling(true,true);
elseif strcmp(site,'prostate') && strcmp(planAlias,'COWC_noAngles')
    sampling = prostateAngularSampling(false,true);
else
    sampling.sampling_scen_mode = 'impScen_permuted5_truncated';
    sampling.sampling_ctActive = strcmp(site,'prostate') || ...
        strcmp(site,'head_and_neck');
    sampling.sampling_setupActive = true;
    sampling.sampling_rangeActive = false;
    sampling.sampling_gantryActive = false;
    sampling.sampling_couchActive = false;
    switch site
        case 'breast'
            sampling.sampling_shiftSD = [4 8 6];
        case 'prostate'
            sampling.sampling_shiftSD = [5 5 5];
        case 'head_and_neck'
            sampling.sampling_shiftSD = [3 3 3];
    end
    sampling.sampling_wcSigma = 1.5;
end

end

function sampling = prostateAngularSampling(withAngles,withRange)

sampling = struct();
sampling.sampling_linkToOptimization = true;
sampling.sampling_scen_mode = 'random';
sampling.sampling_ctActive = false;
sampling.sampling_ctReferenceScenId = 1;
sampling.sampling_setupActive = true;
sampling.sampling_shiftSD = [5 5 5];
sampling.sampling_wcSigma = 1.5;
sampling.sampling_rangeActive = withRange;
if withRange
    sampling.sampling_rangeAbsSD = 1;
    sampling.sampling_rangeRelSD = 3.5;
end
sampling.sampling_gantryActive = withAngles;
sampling.sampling_couchActive = withAngles;
if withAngles
    sampling.sampling_gantryAngleSD = 1;
    sampling.sampling_couchAngleSD = 1;
else
    sampling.sampling_gantryAngleSD = 0;
    sampling.sampling_couchAngleSD = 0;
end

end

function analysis = analysisDefaults()

analysis = struct();
analysis.evaluationMode = 'total';
analysis.gammaWindow = [0 1];
analysis.gammaCriteria = [3 3];
analysis.robustnessCriteria = [5 5];
analysis.robustnessTargetMode = 'include';
analysis.robustnessTargets = {'CTV'};

end

function seed = randomSeedForPlan(site,planAlias)

if strcmp(site,'prostate') && ...
        any(strcmp(planAlias,{'PTV','PTV_noAngles','COWC', ...
        'COWC_noAngles'}))
    seed = 9999;
else
    seed = [];
end

end

function rootPath = userDataRoot()

specFolder = fileparts(mfilename('fullpath'));
sharedFolder = fileparts(specFolder);
macroRoot = fileparts(sharedFolder);
rootPath = fileparts(macroRoot);

end
