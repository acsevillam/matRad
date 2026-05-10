function [pln_interval,dij_intervalContext] = matRad_calcDoseIntervalStreamingCore(ct,cst,stf,pln,cfg,intervalMode)
% matRad_calcDoseIntervalStreamingCore shared streaming implementation for dose interval methods
%
% call
%   [pln_interval,dij_intervalContext] = matRad_calcDoseIntervalStreamingCore(ct,cst,stf,pln,cfg,intervalMode)
%
% input
%   ct:           matRad ct struct
%   cst:          matRad cst cell array
%   stf:          matRad steering information struct
%   pln:          matRad pln struct with robust scenario model
%   cfg:          dose interval streaming configuration struct
%   intervalMode: interval method identifier, either 'INTERVAL2' or 'INTERVAL3'
%
% output
%   pln_interval:        plan struct containing propOpt.dij_interval
%   dij_intervalContext: lightweight single-scenario dij context for
%                        interval fluence optimization
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    cfg = struct();
end
if nargin < 6
    intervalMode = [];
end

matRad_cfg = MatRad_Config.instance();
timer = tic;
cacheContext = [];

try
    [cfg,intervalMode] = normalizeIntervalStreamingConfig(cfg,intervalMode,matRad_cfg);
    scenarioInfo = matRad_resolveStreamingScenarioSelection(ct,pln,cfg,matRad_cfg);
    [provider,dijForResolve] = matRad_initializeScenarioRowProvider( ...
        ct,cst,stf,pln,cfg,scenarioInfo,matRad_cfg,'dose interval');

    ctx = matRad_resolveScenarioDoseInputs(ct,cst,pln,dijForResolve,cfg, ...
        intervalMode,matRad_cfg);
    ctx.scenarioIds = scenarioInfo.scenarioIds(:);
    if isfield(dijForResolve,'doseGrid')
        ctx.doseGrid = dijForResolve.doseGrid;
    end
    if isfield(dijForResolve,'ctGrid')
        ctx.ctGrid = dijForResolve.ctGrid;
    end
    cfg = ctx.cfg;
    cfg.IntervalMode = intervalMode;
    cfg.SecondPassStrategy = provider.secondPassStrategy;
    cfg.KeepCache = provider.keepCache;
    cfg.CacheRoot = provider.cacheRoot;
    quantity = ctx.quantity;

    matRad_cfg.dispInfo(['matRad: Calculating %s dose interval using ', ...
        'streaming scenario rows for quantity ''%s'' and %d scenarios.\n'], ...
        intervalMode,quantity.name,numel(ctx.scenarioDijIx));

    targetBatches = matRad_makeScenarioDoseRowBatches(ctx.targetRows, ...
        matRad_resolveScenarioDoseBatchSize(numel(ctx.targetRows), ...
        numel(ctx.scenarioDijIx),ctx.numBixels,cfg));
    oarBatches = matRad_makeScenarioDoseRowBatches(ctx.oarRows, ...
        matRad_resolveScenarioDoseBatchSize(numel(ctx.oarRows), ...
        numel(ctx.scenarioDijIx),ctx.numBixels,cfg));

    needsTargetSecondPass = ~isempty(ctx.targetRows) && ...
        strcmp(cfg.RadiusMode,'extreme');
    needsOARSecondPass = ~isempty(ctx.oarRows) && strcmp(intervalMode,'INTERVAL3');
    needsDiskCache = strcmp(cfg.SecondPassStrategy,'disk') && ...
        (needsTargetSecondPass || needsOARSecondPass);

    if needsOARSecondPass
        guardIntervalStreamingOARRadiusFactorMemory(numel(ctx.scenarioDijIx), ...
            ctx.numBixels,cfg,matRad_cfg);
    end

    if needsDiskCache
        signatureExtras = struct();
        signatureExtras.targetRows = ctx.targetRows(:);
        signatureExtras.oarRows = ctx.oarRows(:);
        cacheContext = matRad_initializeScenarioDoseCache(cfg,ctx,quantity, ...
            stf,matRad_cfg,'dose interval',signatureExtras);
        provider.cacheContext = cacheContext;
    end

    firstPassWork = struct();
    firstPassWork.primaryBatches = targetBatches;
    firstPassWork.secondaryBatches = oarBatches;
    firstPassWork.cachePrimaryRows = needsTargetSecondPass && ...
        strcmp(cfg.SecondPassStrategy,'disk');
    firstPassWork.cacheSecondaryRows = needsOARSecondPass && ...
        strcmp(cfg.SecondPassStrategy,'disk');

    dij_interval = initializeIntervalStreamingStruct(ctx,quantity,cfg,intervalMode);
    [dij_interval,provider] = accumulateIntervalStreamingFirstPass( ...
        dij_interval,provider,ctx,quantity,cfg,firstPassWork,matRad_cfg);

    if needsTargetSecondPass
        [dij_interval,provider] = accumulateIntervalStreamingTargetExtremeRadius( ...
            dij_interval,provider,ctx,quantity,cfg,targetBatches,matRad_cfg);
    end

    if needsOARSecondPass
        dij_interval = accumulateIntervalStreamingOARRadiusFactors( ...
            dij_interval,provider,ctx,quantity,cfg,oarBatches,matRad_cfg);
    end

    if cfg.CollectTiming
        dij_interval.timing.totalSeconds = toc(timer);
    end

    if needsDiskCache && cfg.KeepCache
        dij_interval.cacheDir = cacheContext.runDir;
    end

    pln_interval = pln;
    if ~isfield(pln_interval,'propOpt') || ~isstruct(pln_interval.propOpt)
        pln_interval.propOpt = struct();
    end
    pln_interval.propOpt.dij_interval = dij_interval;
    dij_intervalContext = buildIntervalStreamingDijContext(dijForResolve, ...
        dij_interval,quantity,cfg);
    pln_interval.multScen = dij_intervalContext.scenarioModel;

    matRad_cleanupScenarioDoseCache(cacheContext,cfg.KeepCache);
    matRad_cfg.dispInfo('matRad: Finished %s streaming dose interval calculation in %.2f s.\n', ...
        intervalMode,toc(timer));
catch ME
    if ~isempty(cacheContext)
        matRad_cleanupScenarioDoseCache(cacheContext,false);
    end
    rethrow(ME);
end

end

function [cfg,intervalMode] = normalizeIntervalStreamingConfig(cfg,intervalMode,matRad_cfg)
if isempty(cfg)
    cfg = struct();
elseif ~isstruct(cfg)
    matRad_cfg.dispError('Dose interval streaming configuration must be a struct.');
end

intervalMode = normalizeIntervalStreamingText(intervalMode,'intervalMode',matRad_cfg);
intervalMode = upper(intervalMode);
if ~any(strcmp(intervalMode,{'INTERVAL2','INTERVAL3'}))
    matRad_cfg.dispError('intervalMode must be ''INTERVAL2'' or ''INTERVAL3''.');
end
cfg.IntervalMode = intervalMode;

cfg = matRad_normalizeScenarioDoseStreamingConfig(cfg,matRad_cfg, ...
    fullfile(matRad_cfg.matRadRoot,'cache','doseInterval'));
end

function dij_interval = initializeIntervalStreamingStruct(ctx,quantity,cfg,intervalMode)
dij_interval = struct();
dij_interval.center = sparse(ctx.numVoxels,ctx.numBixels);
dij_interval.radius = sparse(ctx.numBixels,ctx.numBixels);
dij_interval.targetSubIx = ctx.targetRows(:);
dij_interval.OARSubIx = ctx.oarRows(:);
dij_interval.quantity = quantity.name;
dij_interval.quantityField = quantity.field;
dij_interval.quantityScale = quantity.scale;
dij_interval.optimizationQuantity = quantity.optimizationQuantity;
dij_interval.refScen = cfg.refScen;
dij_interval.scenarioDijIx = ctx.scenarioDijIx(:);
dij_interval.scenarioCtScenIds = ctx.scenarioCtScenIds(:);
dij_interval.scenarioWeights = ctx.scenarioWeights(:);
dij_interval.intervalMode = intervalMode;
dij_interval.radiusMode = cfg.RadiusMode;
dij_interval.precomputeMode = 'streaming';
dij_interval.secondPassStrategy = cfg.SecondPassStrategy;

if cfg.CollectTiming
    dij_interval.timing = matRad_initializeDoseIntervalTiming('interval', ...
        intervalMode,cfg,numel(ctx.targetRows),numel(ctx.oarRows), ...
        numel(ctx.scenarioDijIx),ctx.numBixels);
end

if strcmp(intervalMode,'INTERVAL3')
    dij_interval.OARRadiusRank = zeros(numel(ctx.oarRows),1);
    dij_interval.OARRadiusFactor = cell(numel(ctx.oarRows),1);
end
end

function [dij_interval,provider] = accumulateIntervalStreamingFirstPass( ...
    dij_interval,provider,ctx,quantity,cfg,firstPassWork,matRad_cfg)
numScenarios = numel(ctx.scenarioDijIx);
targetBatches = firstPassWork.primaryBatches;
oarBatches = firstPassWork.secondaryBatches;
cacheTargetRows = firstPassWork.cachePrimaryRows;
cacheOARRows = firstPassWork.cacheSecondaryRows;
targetCenter = sparse(numel(ctx.targetRows),ctx.numBixels);
oarCenter = sparse(numel(ctx.oarRows),ctx.numBixels);
targetSecondMoment = sparse(ctx.numBixels,ctx.numBixels);

matRad_cfg.dispInfo(['matRad: Streaming first pass over %d scenario(s), ', ...
    '%d target voxels, and %d OAR voxels.\n'],numScenarios, ...
    numel(ctx.targetRows),numel(ctx.oarRows));

for s = 1:numScenarios
    matRad_cfg.dispInfo('matRad: Streaming first pass scenario %d/%d.\n', ...
        s,numScenarios);
    [provider,source] = matRad_beginScenarioDoseRowsProvider(provider,ctx,quantity,s,matRad_cfg);
    scenarioWeight = ctx.scenarioWeights(s);

    for b = 1:numel(targetBatches)
        batch = targetBatches{b};
        matRad_logScenarioDoseBatchProgress(matRad_cfg,cfg,'Target center', ...
            b,numel(targetBatches));
        rows = matRad_getScenarioDoseProviderRows(source,ctx.scenarioMaps{s}, ...
            batch.rows,matRad_cfg);
        targetCenter(batch.localIx,:) = targetCenter(batch.localIx,:) + ...
            scenarioWeight .* rows;
        if strcmp(cfg.RadiusMode,'std')
            targetSecondMoment = targetSecondMoment + rows' * ...
                (scenarioWeight .* rows);
        end
        if cacheTargetRows
            matRad_writeScenarioDoseCacheBlock(provider.cacheContext,s, ...
                intervalStreamingTargetCacheKind(),b,batch.rows,rows);
        end
    end

    for b = 1:numel(oarBatches)
        batch = oarBatches{b};
        matRad_logScenarioDoseBatchProgress(matRad_cfg,cfg,'OAR center', ...
            b,numel(oarBatches));
        rows = matRad_getScenarioDoseProviderRows(source,ctx.scenarioMaps{s}, ...
            batch.rows,matRad_cfg);
        oarCenter(batch.localIx,:) = oarCenter(batch.localIx,:) + ...
            scenarioWeight .* rows;
        if cacheOARRows
            matRad_writeScenarioDoseCacheBlock(provider.cacheContext,s, ...
                intervalStreamingOARCacheKind(),b,batch.rows,rows);
        end
    end

    [provider,source] = matRad_endScenarioDoseRowsProvider(provider,source);
end

if ~isempty(ctx.targetRows)
    dij_interval.center(ctx.targetRows,:) = targetCenter;
end
if ~isempty(ctx.oarRows)
    dij_interval.center(ctx.oarRows,:) = oarCenter;
end
if strcmp(cfg.RadiusMode,'std')
    dij_interval.radius = targetSecondMoment;
end
end

function [dij_interval,provider] = accumulateIntervalStreamingTargetExtremeRadius( ...
    dij_interval,provider,ctx,quantity,cfg,targetBatches,matRad_cfg)
targetDeltaRows = sparse(numel(ctx.targetRows),ctx.numBixels);
numScenarios = numel(ctx.scenarioDijIx);

matRad_cfg.dispInfo('matRad: Streaming second pass for extreme target radius.\n');
for s = 1:numScenarios
    source = [];
    if strcmp(cfg.SecondPassStrategy,'recompute')
        [provider,source] = matRad_beginScenarioDoseRowsProvider(provider,ctx,quantity,s,matRad_cfg);
    end

    for b = 1:numel(targetBatches)
        batch = targetBatches{b};
        matRad_logScenarioDoseBatchProgress(matRad_cfg,cfg,'Target extreme radius', ...
            b,numel(targetBatches));
        if strcmp(cfg.SecondPassStrategy,'disk')
            rows = matRad_readScenarioDoseCacheBlock(provider.cacheContext,s, ...
                intervalStreamingTargetCacheKind(),b,batch.rows);
        else
            rows = matRad_getScenarioDoseProviderRows(source,ctx.scenarioMaps{s}, ...
                batch.rows,matRad_cfg);
        end
        centerRows = dij_interval.center(batch.rows,:);
        targetDeltaRows(batch.localIx,:) = max(targetDeltaRows(batch.localIx,:), ...
            abs(rows - centerRows));
    end

    if strcmp(cfg.SecondPassStrategy,'recompute')
        [provider,source] = matRad_endScenarioDoseRowsProvider(provider,source);
    end
end

deltaSquared = full(sum(targetDeltaRows.^2,1));
centerRows = dij_interval.center(ctx.targetRows,:);
dij_interval.radius = centerRows' * centerRows + ...
    spdiags(deltaSquared(:),0,ctx.numBixels,ctx.numBixels);
end

function [dij_interval,provider] = accumulateIntervalStreamingOARRadiusFactors( ...
    dij_interval,provider,ctx,quantity,cfg,oarBatches,matRad_cfg)
numScenarios = numel(ctx.scenarioDijIx);
intervalOffset = 0;

matRad_cfg.dispInfo('matRad: Streaming second pass for INTERVAL3 OAR radius factors.\n');
for b = 1:numel(oarBatches)
    batch = oarBatches{b};
    scenarioRows = cell(numScenarios,1);
    matRad_logScenarioDoseBatchProgress(matRad_cfg,cfg,'OAR radius factor', ...
        b,numel(oarBatches));

    for s = 1:numScenarios
        if strcmp(cfg.SecondPassStrategy,'disk')
            scenarioRows{s} = matRad_readScenarioDoseCacheBlock(provider.cacheContext, ...
                s,intervalStreamingOARCacheKind(),b,batch.rows);
        else
            [provider,source] = matRad_beginScenarioDoseRowsProvider(provider,ctx,quantity,s,matRad_cfg);
            scenarioRows{s} = matRad_getScenarioDoseProviderRows(source,ctx.scenarioMaps{s}, ...
                batch.rows,matRad_cfg);
            [provider,source] = matRad_endScenarioDoseRowsProvider(provider,source);
        end
    end

    for localIx = 1:numel(batch.rows)
        centeredMatrixRows = cell(numScenarios,1);
        centerRow = dij_interval.center(batch.rows(localIx),:);
        for s = 1:numScenarios
            centeredMatrixRows{s} = scenarioRows{s}(localIx,:) - centerRow;
        end
        centeredScenarioMatrix = vertcat(centeredMatrixRows{:});
        [factor,rank] = buildIntervalStreamingRadiusFactor( ...
            centeredScenarioMatrix,ctx.scenarioWeights,cfg,ctx.numBixels);
        dij_interval.OARRadiusRank(intervalOffset + localIx) = rank;
        dij_interval.OARRadiusFactor{intervalOffset + localIx} = factor;
    end

    intervalOffset = intervalOffset + numel(batch.rows);
end
end

function dij_intervalContext = buildIntervalStreamingDijContext(dij,dij_interval,quantity,cfg)
metadataFields = {'doseGrid','ctGrid','totalNumOfBixels','numOfBeams', ...
    'beamNum','rayNum','bixelNum','numParticlesPerMU','minMU','maxMU', ...
    'RBE','RBE_models','ax','bx','doseWeightingThreshold','machine','meta'};
dij_intervalContext = struct();
for f = 1:numel(metadataFields)
    fieldName = metadataFields{f};
    if isfield(dij,fieldName)
        dij_intervalContext.(fieldName) = dij.(fieldName);
    end
end

numBixels = size(dij_interval.center,2);
dij_intervalContext.totalNumOfBixels = numBixels;
if ~isfield(dij_intervalContext,'beamNum') || ...
        numel(dij_intervalContext.beamNum) ~= numBixels
    dij_intervalContext.beamNum = ones(numBixels,1);
else
    dij_intervalContext.beamNum = dij_intervalContext.beamNum(:);
end
if ~isfield(dij_intervalContext,'numOfBeams') || ...
        isempty(dij_intervalContext.numOfBeams)
    dij_intervalContext.numOfBeams = max(dij_intervalContext.beamNum);
end

dij_intervalContext.physicalDose = cell(1,1,1);
dij_intervalContext.physicalDose{1} = dij_interval.center;
if ~strcmp(quantity.field,'physicalDose')
    dij_intervalContext.(quantity.field) = cell(1,1,1);
    dij_intervalContext.(quantity.field){1} = dij_interval.center;
end
dij_intervalContext.numOfScenarios = 1;
dij_intervalContext.scenarioModel = matRad_buildStreamingOptimizationScenarioModel( ...
    cfg.refScen);
dij_intervalContext.intervalQuantity = quantity.name;
dij_intervalContext.intervalQuantityField = quantity.field;
end

function [factor,rank] = buildIntervalStreamingRadiusFactor(centeredScenarioMatrix, ...
    scenarioWeights,cfg,numBixels)
if strcmp(cfg.RadiusMode,'extreme')
    [factor,rank] = buildIntervalStreamingExtremeRadiusFactor( ...
        centeredScenarioMatrix,scenarioWeights,numBixels);
else
    [factor,rank] = buildIntervalStreamingStdRadiusFactor(centeredScenarioMatrix, ...
        scenarioWeights,cfg,numBixels);
end
end

function [factor,rank] = buildIntervalStreamingStdRadiusFactor(centeredScenarioMatrix, ...
    scenarioWeights,cfg,numBixels)
tolerance = 1e-12;

kMax = min(cfg.KMax,numBixels);
if kMax == 0
    [factor,rank] = zeroRankIntervalStreamingRadiusFactor(numBixels);
    return;
end

positiveWeightIx = scenarioWeights(:) > 0;
scenarioWeights = scenarioWeights(positiveWeightIx);
centeredScenarioMatrix = centeredScenarioMatrix(positiveWeightIx,:);

activeBixels = find(any(centeredScenarioMatrix,1));
if isempty(activeBixels)
    [factor,rank] = zeroRankIntervalStreamingRadiusFactor(numBixels);
    return;
end

centeredScenarioMatrix = centeredScenarioMatrix(:,activeBixels);
numScenarios = size(centeredScenarioMatrix,1);
weightedScenarioMatrix = spdiags(sqrt(scenarioWeights),0, ...
    numScenarios,numScenarios) * centeredScenarioMatrix;

if nnz(weightedScenarioMatrix) == 0 || ...
        norm(weightedScenarioMatrix,'fro') <= tolerance
    [factor,rank] = zeroRankIntervalStreamingRadiusFactor(numBixels);
    return;
end

scenarioGram = weightedScenarioMatrix * weightedScenarioMatrix';
scenarioGram = full(0.5 .* (scenarioGram + scenarioGram'));
[leftEigenVectors,eigenValueMatrix] = eig(scenarioGram);
eigenValues = real(diag(eigenValueMatrix));
eigenValues(abs(eigenValues) <= tolerance) = 0;
eigenValues(eigenValues < 0) = 0;
[eigenValues,sortIx] = sort(eigenValues,'descend');
leftEigenVectors = leftEigenVectors(:,sortIx);

positiveIx = find(eigenValues > tolerance);
if isempty(positiveIx)
    [factor,rank] = zeroRankIntervalStreamingRadiusFactor(numBixels);
    return;
end

eigenValues = eigenValues(positiveIx);
leftEigenVectors = leftEigenVectors(:,positiveIx);

if strcmp(cfg.KMode,'dynamic')
    rank = selectIntervalStreamingFactorRank(eigenValues,kMax, ...
        cfg.RetentionThreshold,tolerance);
else
    rank = min(kMax,numel(eigenValues));
end

if rank == 0
    [factor,rank] = zeroRankIntervalStreamingRadiusFactor(numBixels);
    return;
end

eigenValues = eigenValues(1:rank);
leftEigenVectors = leftEigenVectors(:,1:rank);
rightFactorBasis = weightedScenarioMatrix' * leftEigenVectors;
rightFactorBasis = bsxfun(@rdivide,rightFactorBasis, ...
    sqrt(eigenValues(:))');
factorActive = bsxfun(@times,rightFactorBasis,sqrt(eigenValues(:))');

factor = sparse(double(numBixels),double(rank));
factor(activeBixels,:) = sparse(factorActive);
end

function [factor,rank] = buildIntervalStreamingExtremeRadiusFactor(centeredScenarioMatrix, ...
    scenarioWeights,numBixels)
tolerance = 1e-12;

positiveWeightIx = scenarioWeights(:) > 0;
centeredScenarioMatrix = centeredScenarioMatrix(positiveWeightIx,:);
if isempty(centeredScenarioMatrix)
    [factor,rank] = zeroRankIntervalStreamingRadiusFactor(numBixels);
    return;
end

delta = max(abs(centeredScenarioMatrix),[],1);
delta(abs(delta) <= tolerance) = 0;
activeBixels = find(delta);

rank = numel(activeBixels);
if rank == 0
    [factor,rank] = zeroRankIntervalStreamingRadiusFactor(numBixels);
    return;
end

deltaValues = full(delta(activeBixels(:)));
factor = sparse(activeBixels(:),(1:rank)',deltaValues(:), ...
    double(numBixels),double(rank));
end

function rank = selectIntervalStreamingFactorRank(eigenValues,kMax,retentionThreshold,tolerance)
maxRank = min(kMax,numel(eigenValues));
if maxRank == 0
    rank = 0;
    return;
end

rankEigenValues = eigenValues(1:maxRank);
energyScale = max(abs(rankEigenValues));
if ~isfinite(energyScale) || energyScale <= tolerance
    rank = maxRank;
    return;
end

relativeEnergy = (rankEigenValues ./ energyScale).^2;
totalEnergy = sum(relativeEnergy);
if ~isfinite(totalEnergy) || totalEnergy <= 0
    rank = maxRank;
    return;
end

thresholdEnergy = retentionThreshold * totalEnergy;
rank = find(cumsum(relativeEnergy) >= thresholdEnergy,1,'first');
if isempty(rank)
    rank = maxRank;
end
end

function [factor,rank] = zeroRankIntervalStreamingRadiusFactor(numBixels)
rank = 0;
factor = sparse(numBixels,0);
end

function guardIntervalStreamingOARRadiusFactorMemory(numScenarios,numBixels,cfg,matRad_cfg)
memoryLimitMB = matRad_resolveScenarioDoseMemoryLimitMB(cfg);
estimatedMB = estimateIntervalStreamingOARRadiusFactorMemoryMB(numScenarios,numBixels,cfg);

matRad_cfg.dispInfo(['matRad: Estimated INTERVAL3 OAR radius factor memory ', ...
    'per voxel is %.2f MB (limit %.2f MB).\n'],estimatedMB,memoryLimitMB);

if estimatedMB > memoryLimitMB
    matRad_cfg.dispError(['INTERVAL3 OAR radius factor estimated memory per ', ...
        'voxel is %.2f MB, which exceeds MemoryLimitMB %.2f MB. Increase ', ...
        'cfg.MemoryLimitMB, reduce cfg.KMax, reduce the number of scenarios ', ...
        'or bixels, or use INTERVAL2.'],estimatedMB,memoryLimitMB);
end
end

function estimatedMB = estimateIntervalStreamingOARRadiusFactorMemoryMB(numScenarios,numBixels,cfg)
bytesPerDouble = 8;
kMax = min([cfg.KMax,numBixels,numScenarios]);

scenarioMatrixBytes = 3 * numScenarios * numBixels * bytesPerDouble;
scenarioGramBytes = 3 * numScenarios^2 * bytesPerDouble;
factorBytes = (3*numBixels*kMax + kMax^2) * bytesPerDouble;

estimatedMB = (scenarioMatrixBytes + scenarioGramBytes + ...
    factorBytes) / 1e6;
end

function textValue = normalizeIntervalStreamingText(value,name,matRad_cfg)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value) || isempty(value)
    matRad_cfg.dispError('%s must be a non-empty character vector.',name);
end
textValue = value;
end

function kind = intervalStreamingTargetCacheKind()
kind = 'target';
end

function kind = intervalStreamingOARCacheKind()
kind = 'oar';
end
