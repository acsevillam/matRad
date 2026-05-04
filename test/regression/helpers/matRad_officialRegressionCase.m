function result = matRad_officialRegressionCase(caseName)
% Deterministic regression cases shared by local and official checkouts.

switch caseName
    case 'scenario_generation'
        result = scenarioGenerationCase();
    case 'dose_engine_initialization'
        result = doseEngineInitializationCase();
    case 'pencil_beam_dose'
        result = pencilBeamDoseCase();
    case 'particle_bio_pencil_beam_dose'
        result = particleBioPencilBeamDoseCase();
    otherwise
        error('matRad:RegressionCaseFailed', ...
            'Unknown regression case "%s".',caseName);
end

end

function result = scenarioGenerationCase()
    ct.numOfCtScen = 3;
    ctScenProb = [1 0.2; 2 0.3; 3 0.5];

    result.nominal = captureScenario(configureCommonScenario( ...
        createScenarioModel('nomScen',ct),ctScenProb), ...
        'nominal_three_ct');

    wcModel = configureCommonScenario(createScenarioModel('wcScen',ct),ctScenProb);
    result.worstCaseCombined = captureScenario(wcModel, ...
        'worst_case_combined_range');

    wcShiftModel = configureCommonScenario(createScenarioModel('wcScen',ct),ctScenProb);
    setIfProp(wcShiftModel,'combineRange',false);
    setIfProp(wcShiftModel,'combinations','shift');
    result.worstCaseShift = captureScenario(wcShiftModel, ...
        'worst_case_shift_combinations');

    importanceModel = configureCommonScenario(createScenarioModel('impScen',ct),ctScenProb);
    setIfProp(importanceModel,'combinations','none');
    setIfProp(importanceModel,'combineRange',true);
    result.importance = captureScenario(importanceModel, ...
        'importance_small_grid');

    randomModel = configureCommonScenario(createScenarioModel('rndScen',ct),ctScenProb);
    setIfProp(randomModel,'includeNominalScenario',true);
    setRegressionSeed(7193);
    setIfProp(randomModel,'nSamples',5);
    result.random = captureScenario(randomModel, ...
        'random_seeded_samples');
end

function result = doseEngineInitializationCase()
    [ct,cst,stf,pln] = doseEngineFixture();
    engine = DoseEngines.matRad_RegressionProbeEngine(pln);
    dij = engine.calcDoseInfluence(ct,cst,stf);

    result.dij = captureStructFields(dij,{ ...
        'ctGrid','doseGrid','numOfBeams','numOfScenarios', ...
        'numOfRaysPerBeam','totalNumOfBixels','totalNumOfRays', ...
        'bixelNum','rayNum','beamNum','minMU','maxMU', ...
        'numOfParticlesPerMU','regressionProbe'});
    result.scenario = captureScenario(pln.multScen,'engine_mult_scen');
end

function result = pencilBeamDoseCase()
    result.photons = pencilBeamDoseForTestData('photons_testData.mat','SVDPB');
    result.protons = pencilBeamDoseForTestData('protons_testData.mat','HongPB');
    result.helium = pencilBeamDoseForTestData('helium_testData.mat','HongPB');
    result.carbon = pencilBeamDoseForTestData('carbon_testData.mat','HongPB');
end

function result = particleBioPencilBeamDoseCase()
    result.protonsConstRBE = pencilBeamDoseForTestData( ...
        'protons_testData.mat','HongPB','RBExD','constRBE');
    result.protonsMCN = pencilBeamDoseForTestData( ...
        'protons_testData.mat','HongPB','RBExD','MCN');
    result.protonsWED = pencilBeamDoseForTestData( ...
        'protons_testData.mat','HongPB','RBExD','WED');
    result.heliumHEL = pencilBeamDoseForTestData( ...
        'helium_testData.mat','HongPB','RBExD','HEL');
    result.carbonLEM = pencilBeamDoseForTestData( ...
        'carbon_testData.mat','HongPB','RBExD','LEM');
end

function result = pencilBeamDoseForTestData(fileName,engineName,quantityOpt,modelName)
    if nargin < 3
        quantityOpt = 'physicalDose';
    end
    if nargin < 4
        modelName = 'none';
    end

    data = load(fullfile(pwd,'test','testData',fileName),'ct','cst','pln','stf');
    ct = data.ct;
    cst = data.cst;
    pln = data.pln;
    stf = data.stf;

    pln.propDoseCalc.engine = engineName;
    pln.multScen = createScenarioModel('nomScen',ct);
    pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

    dij = matRad_calcDoseInfluence(ct,cst,stf,pln);
    assertDoseDijIsPhysicallyPlausible(dij);

    weights = deterministicWeights(dij.totalNumOfBixels);

    result.fileName = fileName;
    result.radiationMode = pln.radiationMode;
    result.engineName = engineName;
    result.quantityOpt = quantityOpt;
    result.modelName = modelName;
    result.bioParam = captureObjectProperties(pln.bioParam,{ ...
        'radiationMode','quantityOpt','quantityVis','model','RBE','bioOpt'});
    result.dij = captureStructFields(dij,{ ...
        'ctGrid','doseGrid','numOfBeams','numOfScenarios', ...
        'numOfRaysPerBeam','totalNumOfBixels','totalNumOfRays', ...
        'bixelNum','rayNum','beamNum','physicalDose','mLETDose', ...
        'mAlphaDose','mSqrtBetaDose','RBE','ax','bx','vTissueIndex'});
    result.weights = weights;
    result.forward = calculateForwardDoseSummary(dij,weights);
end

function weights = deterministicWeights(numOfBixels)
    weights = (1:numOfBixels)';
    weights = 1 + weights ./ max(1,numOfBixels);
end

function assertDoseDijIsPhysicallyPlausible(dij)
    if ~isfield(dij,'physicalDose') || ~iscell(dij.physicalDose)
        error('matRad:RegressionCaseFailed', ...
            'Pencil-beam dij has no physicalDose cell array.');
    end

    populatedScenarios = find(~cellfun(@isempty,dij.physicalDose(:)));
    if isempty(populatedScenarios)
        error('matRad:RegressionCaseFailed', ...
            'Pencil-beam dij has no populated physicalDose scenarios.');
    end

    for i = 1:numel(populatedScenarios)
        scenarioIx = populatedScenarios(i);
        doseMatrix = full(dij.physicalDose{scenarioIx});
        if size(doseMatrix,1) ~= dij.doseGrid.numOfVoxels || ...
                size(doseMatrix,2) ~= dij.totalNumOfBixels
            error('matRad:RegressionCaseFailed', ...
                'physicalDose{%d} has inconsistent dimensions.',scenarioIx);
        end
        if any(~isfinite(doseMatrix(:)))
            error('matRad:RegressionCaseFailed', ...
                'physicalDose{%d} contains non-finite values.',scenarioIx);
        end
        if any(doseMatrix(:) < -1e-12)
            error('matRad:RegressionCaseFailed', ...
                'physicalDose{%d} contains negative dose values.',scenarioIx);
        end
        if ~any(doseMatrix(:) > 0)
            error('matRad:RegressionCaseFailed', ...
                'physicalDose{%d} is identically zero.',scenarioIx);
        end
    end

    if numel(dij.beamNum) ~= dij.totalNumOfBixels || ...
            numel(dij.rayNum) ~= dij.totalNumOfBixels || ...
            numel(dij.bixelNum) ~= dij.totalNumOfBixels
        error('matRad:RegressionCaseFailed', ...
            'Dij beam/ray/bixel bookkeeping length is inconsistent.');
    end
end

function forward = calculateForwardDoseSummary(dij,weights)
    fullScenIx = 1;
    physicalDose = reshape(full(dij.physicalDose{fullScenIx} * weights), ...
        dij.doseGrid.dimensions);

    forward.physicalDose = normalizeValue(physicalDose);
    forward.physicalDoseSum = sum(physicalDose(:));
    forward.physicalDoseMax = max(physicalDose(:));
    forward.physicalDoseMin = min(physicalDose(:));

    if isfield(dij,'RBE') && isscalar(dij.RBE)
        forward.RBExD = physicalDose .* dij.RBE;
        forward.RBExDSum = sum(forward.RBExD(:));
    end

    if isfield(dij,'mLETDose') && ~isempty(dij.mLETDose{fullScenIx})
        mLETDose = reshape(full(dij.mLETDose{fullScenIx} * weights), ...
            dij.doseGrid.dimensions);
        forward.mLETDose = normalizeValue(mLETDose);
        forward.mLETDoseSum = sum(mLETDose(:));
    end

    if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose') && ...
            ~isempty(dij.mAlphaDose{fullScenIx}) && ...
            ~isempty(dij.mSqrtBetaDose{fullScenIx})
        mAlphaDose = reshape(full(dij.mAlphaDose{fullScenIx} * weights), ...
            dij.doseGrid.dimensions);
        mSqrtBetaDose = reshape(full(dij.mSqrtBetaDose{fullScenIx} * weights), ...
            dij.doseGrid.dimensions);
        effect = mAlphaDose + mSqrtBetaDose.^2;
        forward.mAlphaDose = normalizeValue(mAlphaDose);
        forward.mSqrtBetaDose = normalizeValue(mSqrtBetaDose);
        forward.effect = normalizeValue(effect);
        forward.effectSum = sum(effect(:));

        if isfield(dij,'ax') && isfield(dij,'bx')
            ctScenId = 1;
            ax = dij.ax{ctScenId};
            bx = dij.bx{ctScenId};
            rbexd = zeros(size(physicalDose));
            rbexdIx = bx ~= 0 & physicalDose(:) > 0;
            rbexd(rbexdIx) = (sqrt(ax(rbexdIx).^2 + ...
                4 .* bx(rbexdIx) .* effect(rbexdIx)) - ax(rbexdIx)) ./ ...
                (2 .* bx(rbexdIx));
            forward.RBExD = normalizeValue(rbexd);
            forward.RBExDSum = sum(rbexd(:));
        end
    end

    forward = orderfields(forward);
end

function [ct,cst,stf,pln] = doseEngineFixture()
    ct.numOfCtScen = 2;
    ct.cubeDim = [2 3 2];
    ct.resolution = struct('x',1,'y',1,'z',2);
    ct.refScen = 1;

    cst = cell(2,6);
    cst = addStructure(cst,1,'PTV','TARGET',{[1;2;3;4;5],[2;3;4;5;6]});
    cst = addStructure(cst,2,'OAR','OAR',{[7;8;9;10],[8;9;10;11]});

    stf = repmat(struct(),1,2);
    stf(1).machine = 'Generic';
    stf(1).radiationMode = 'photons';
    stf(1).numOfRays = 2;
    stf(1).numOfBixelsPerRay = [1 2];
    stf(1).totalNumOfBixels = 3;
    stf(2).machine = 'Generic';
    stf(2).radiationMode = 'photons';
    stf(2).numOfRays = 1;
    stf(2).numOfBixelsPerRay = 2;
    stf(2).totalNumOfBixels = 2;

    pln.radiationMode = 'photons';
    pln.machine = 'Generic';
    pln.propDoseCalc.doseGrid.resolution = struct('x',1,'y',1,'z',2);
    pln.propDoseCalc.selectVoxelsInScenarios = 'all';

    ctScenProb = [1 0.4; 2 0.6];
    pln.multScen = configureCommonScenario(createScenarioModel('wcScen',ct),ctScenProb);
    setIfProp(pln.multScen,'combinations','shift');
end

function cst = addStructure(cst,rowIx,name,type,voxelsByCtScenario)
    cst{rowIx,1} = rowIx;
    cst{rowIx,2} = name;
    cst{rowIx,3} = type;
    cst{rowIx,4} = voxelsByCtScenario;
    cst{rowIx,5} = struct('Priority',rowIx);
    cst{rowIx,6} = {};
end

function model = createScenarioModel(modelName,ct)
    if exist('matRad_createScenarioModel','file') == 2
        model = matRad_createScenarioModel(ct,modelName);
        return;
    end

    switch modelName
        case 'nomScen'
            model = matRad_NominalScenario(ct);
        case 'wcScen'
            model = matRad_WorstCaseScenarios(ct);
        case 'impScen'
            model = matRad_ImportanceScenarios(ct);
        case 'rndScen'
            model = matRad_RandomScenarios(ct);
        otherwise
            error('matRad:RegressionCaseFailed', ...
                'Cannot create scenario model "%s".',modelName);
    end
end

function model = configureCommonScenario(model,ctScenProb)
    setIfProp(model,'ctScenProb',ctScenProb);
    setIfProp(model,'shiftSD',[1.25 2.5 3.75]);
    setIfProp(model,'rangeAbsSD',1.5);
    setIfProp(model,'rangeRelSD',2.75);
    setIfProp(model,'wcSigma',1.2);
end

function setIfProp(obj,propName,value)
    if isprop(obj,propName)
        obj.(propName) = value;
    end
end

function data = captureScenario(model,label)
    data.label = label;
    data.className = class(model);
    data.common = captureObjectProperties(model,{ ...
        'name','ctScenProb','numOfCtScen','numOfAvailableCtScen', ...
        'ctScenIx','isoShift','relRangeShift','absRangeShift', ...
        'maxAbsRangeShift','maxRelRangeShift','totNumShiftScen', ...
        'totNumRangeScen','totNumScen','scenForProb','scenProb', ...
        'scenWeight','scenMask','linearMask'});

end

function values = captureObjectProperties(obj,fieldNames)
    values = struct();
    for i = 1:numel(fieldNames)
        fieldName = fieldNames{i};
        if isprop(obj,fieldName)
            values.(fieldName) = normalizeValue(obj.(fieldName));
        end
    end
    values = orderfields(values);
end

function values = captureStructFields(inputStruct,fieldNames)
    values = struct();
    for i = 1:numel(fieldNames)
        fieldName = fieldNames{i};
        if isfield(inputStruct,fieldName)
            values.(fieldName) = normalizeValue(inputStruct.(fieldName));
        end
    end
    values = orderfields(values);
end

function value = normalizeValue(value)
    if isstruct(value)
        for i = 1:numel(value)
            fieldNames = fieldnames(value(i));
            for j = 1:numel(fieldNames)
                value(i).(fieldNames{j}) = normalizeValue(value(i).(fieldNames{j}));
            end
        end
        value = orderfields(value);
    elseif iscell(value)
        for i = 1:numel(value)
            value{i} = normalizeValue(value{i});
        end
    elseif issparse(value)
        value = full(value);
    elseif isstring(value)
        value = cellstr(value);
    end
end

function setRegressionSeed(seed)
    if exist('rng','file') == 2 || exist('rng','builtin') == 5
        rng(seed,'twister');
    else
        rand('seed',seed);
        randn('seed',seed);
    end
end
