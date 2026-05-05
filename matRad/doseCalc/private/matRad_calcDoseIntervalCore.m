function [pln_interval,dij_intervalContext] = matRad_calcDoseIntervalCore(ct,cst,stf,pln,dij,cfg,intervalMode)
% matRad_calcDoseIntervalCore shared implementation for dose interval methods
%
% call
%   [pln_interval,dij_intervalContext] = matRad_calcDoseIntervalCore(ct,cst,stf,pln,dij,cfg,intervalMode)
%
% input
%   ct:           matRad ct struct; multi-CT inputs require pull DVFs when
%                 non-reference CT scenarios are mapped to cfg.refScen
%   cst:          matRad cst cell array
%   stf:          matRad steering information struct
%   pln:          matRad pln struct with a matRad_ScenarioModel in
%                 pln.multScen; optional pln.propOpt.scen4D selects CT
%                 scenarios used for the interval calculation (default: 1)
%   dij:          robust dose influence struct; scenario cells are addressed
%                 by DIJ linear scenario indices from pln.multScen
%   cfg:          dose interval configuration struct
%   intervalMode: interval method identifier, either 'INTERVAL2' or 'INTERVAL3'
%
% output
%   pln_interval:        plan struct containing propOpt.dij_interval
%   dij_intervalContext: lightweight single-scenario dij context for
%                        interval fluence optimization
%
% References
%   -
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

matRad_cfg = MatRad_Config.instance();
timer = tic;

ctx = matRad_resolveDoseIntervalInputs(ct,cst,pln,dij,cfg, ...
    intervalMode,matRad_cfg);
cfg = ctx.cfg;
quantity = ctx.quantity;
scenarioDijIx = ctx.scenarioDijIx;
scenarioCtScenIds = ctx.scenarioCtScenIds;
scenarioWeights = ctx.scenarioWeights;
targetRows = ctx.targetRows;
oarRows = ctx.oarRows;
numVoxels = ctx.numVoxels;
numBixels = ctx.numBixels;
scenarioMaps = ctx.scenarioMaps;

matRad_cfg.dispInfo(['matRad: Calculating %s dose interval for quantity ', ...
    '''%s'' using dij.%s and %d scenarios.\n'],intervalMode, ...
    quantity.name,quantity.field,numel(scenarioDijIx));
matRad_cfg.dispInfo(['matRad: Dose interval reference CT scenario %d, ', ...
    '%d target voxels, %d OAR voxels, %d bixels.\n'],cfg.refScen, ...
    numel(targetRows),numel(oarRows),numBixels);
if quantity.scale ~= 1
    matRad_cfg.dispInfo('matRad: Dose interval quantity scale is %.6g.\n', ...
        quantity.scale);
end

dij_interval = initializeIntervalStruct(numVoxels,numBixels,targetRows,oarRows, ...
    quantity,cfg,scenarioDijIx,scenarioCtScenIds,scenarioWeights,intervalMode);

if ~isempty(targetRows)
    dij_interval = accumulateTargetInterval(dij_interval,quantity,scenarioDijIx, ...
        scenarioWeights,scenarioMaps,targetRows,cfg,numBixels,matRad_cfg);
else
    matRad_cfg.dispInfo('matRad: No target voxels selected for interval target term.\n');
end

if ~isempty(oarRows)
    dij_interval = accumulateOARCenter(dij_interval,quantity,scenarioDijIx, ...
        scenarioWeights,scenarioMaps,oarRows,cfg,numBixels,matRad_cfg);
else
    matRad_cfg.dispInfo('matRad: No OAR voxels selected for interval OAR term.\n');
end

if strcmp(intervalMode,'INTERVAL3') && ~isempty(oarRows)
    guardOARSvdMemory(numel(scenarioDijIx),numBixels,cfg,matRad_cfg);
    dij_interval = accumulateOARSvd(dij_interval,quantity,scenarioDijIx, ...
        scenarioWeights,scenarioMaps,oarRows,cfg,numBixels,matRad_cfg);
end

pln_interval = pln;
if ~isfield(pln_interval,'propOpt') || ~isstruct(pln_interval.propOpt)
    pln_interval.propOpt = struct();
end
pln_interval.propOpt.dij_interval = dij_interval;
dij_intervalContext = buildIntervalDijContext(dij,dij_interval,quantity,cfg);
pln_interval.multScen = dij_intervalContext.scenarioModel;

matRad_cfg.dispInfo('matRad: Finished %s dose interval calculation in %.2f s.\n', ...
    intervalMode,toc(timer));

end

function dij_intervalContext = buildIntervalDijContext(dij,dij_interval,quantity,cfg)
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
dij_intervalContext.scenarioModel = buildIntervalOptimizationScenarioModel( ...
    cfg.refScen);
dij_intervalContext.intervalQuantity = quantity.name;
dij_intervalContext.intervalQuantityField = quantity.field;
end

function scenarioModel = buildIntervalOptimizationScenarioModel(refScen)
components = matRad_createScenarioComponents([0 0 0],0,0,{'ct'});
scenarioValues = zeros(1,numel(components));
ctScenIds = refScen;
scenProb = 1;
scenWeight = 1;
scenForProb = [ctScenIds scenarioValues];
linearMask = [1 1 1];
scenMask = true(1,1,1);

scenarioModel = matRad_NominalScenario();
scenarioModel.setScenarioRealizations(components,scenarioValues,ctScenIds, ...
    scenProb,scenWeight,scenForProb,linearMask,scenMask);
end

function dij_interval = initializeIntervalStruct(numVoxels,numBixels,targetRows,oarRows, ...
    quantity,cfg,scenarioDijIx,scenarioCtScenIds,scenarioWeights,intervalMode)
dij_interval = struct();
dij_interval.center = sparse(numVoxels,numBixels);
dij_interval.radius = sparse(numBixels,numBixels);
dij_interval.targetSubIx = targetRows(:);
dij_interval.OARSubIx = oarRows(:);
dij_interval.quantity = quantity.name;
dij_interval.quantityField = quantity.field;
dij_interval.quantityScale = quantity.scale;
dij_interval.optimizationQuantity = quantity.optimizationQuantity;
dij_interval.refScen = cfg.refScen;
dij_interval.scenarioDijIx = scenarioDijIx(:);
dij_interval.scenarioCtScenIds = scenarioCtScenIds(:);
dij_interval.scenarioWeights = scenarioWeights(:);
dij_interval.intervalMode = intervalMode;

if strcmp(intervalMode,'INTERVAL3')
    dij_interval.k = zeros(numel(oarRows),1);
    dij_interval.U = cell(numel(oarRows),1);
    dij_interval.S = cell(numel(oarRows),1);
    dij_interval.V = cell(numel(oarRows),1);
end
end

function dij_interval = accumulateTargetInterval(dij_interval,quantity,scenarioDijIx, ...
    scenarioWeights,scenarioMaps,targetRows,cfg,numBixels,matRad_cfg)
batchSize = resolveBatchSize(numel(targetRows),numel(scenarioDijIx),numBixels,cfg);
batches = makeBatches(targetRows,batchSize);
matRad_cfg.dispInfo(['matRad: Accumulating target interval center/radius ', ...
    'for %d voxels in %d batches of up to %d voxels.\n'], ...
    numel(targetRows),numel(batches),batchSize);

for b = 1:numel(batches)
    batchTimer = tic;
    rows = batches{b};
    logBatchProgress(matRad_cfg,cfg,'Target interval',b,numel(batches));
    logBatchStart(matRad_cfg,cfg,'Target interval',b,numel(batches),numel(rows));
    centerBlock = sparse(numel(rows),numBixels);
    radiusBlock = sparse(numBixels,numBixels);

    for s = 1:numel(scenarioDijIx)
        logScenarioProgress(matRad_cfg,cfg,'Target interval',b,numel(batches), ...
            s,numel(scenarioDijIx),scenarioDijIx(s),scenarioWeights(s));
        scenarioRows = matRad_getDoseIntervalScenarioRows(quantity,scenarioDijIx(s), ...
            scenarioMaps{s},rows,matRad_cfg);
        centerBlock = centerBlock + scenarioWeights(s) .* scenarioRows;
        radiusBlock = radiusBlock + scenarioRows' * (scenarioWeights(s) .* scenarioRows);
    end

    dij_interval.center(rows,:) = centerBlock;
    dij_interval.radius = dij_interval.radius + radiusBlock;
    logBatchEnd(matRad_cfg,cfg,'Target interval',b,numel(batches),toc(batchTimer));
end
matRad_cfg.dispInfo('matRad: Target interval accumulation finished.\n');
end

function dij_interval = accumulateOARCenter(dij_interval,quantity,scenarioDijIx, ...
    scenarioWeights,scenarioMaps,oarRows,cfg,numBixels,matRad_cfg)
batchSize = resolveBatchSize(numel(oarRows),numel(scenarioDijIx),numBixels,cfg);
batches = makeBatches(oarRows,batchSize);
matRad_cfg.dispInfo(['matRad: Accumulating OAR interval center for %d voxels ', ...
    'in %d batches of up to %d voxels.\n'], ...
    numel(oarRows),numel(batches),batchSize);

for b = 1:numel(batches)
    batchTimer = tic;
    rows = batches{b};
    logBatchProgress(matRad_cfg,cfg,'OAR center',b,numel(batches));
    logBatchStart(matRad_cfg,cfg,'OAR center',b,numel(batches),numel(rows));
    centerBlock = sparse(numel(rows),numBixels);
    for s = 1:numel(scenarioDijIx)
        logScenarioProgress(matRad_cfg,cfg,'OAR center',b,numel(batches), ...
            s,numel(scenarioDijIx),scenarioDijIx(s),scenarioWeights(s));
        scenarioRows = matRad_getDoseIntervalScenarioRows(quantity,scenarioDijIx(s), ...
            scenarioMaps{s},rows,matRad_cfg);
        centerBlock = centerBlock + scenarioWeights(s) .* scenarioRows;
    end
    dij_interval.center(rows,:) = centerBlock;
    logBatchEnd(matRad_cfg,cfg,'OAR center',b,numel(batches),toc(batchTimer));
end
matRad_cfg.dispInfo('matRad: OAR interval center accumulation finished.\n');
end

function dij_interval = accumulateOARSvd(dij_interval,quantity,scenarioDijIx, ...
    scenarioWeights,scenarioMaps,oarRows,cfg,numBixels,matRad_cfg)
batchSize = resolveBatchSize(numel(oarRows),numel(scenarioDijIx),numBixels,cfg);
batches = makeBatches(oarRows,batchSize);
matRad_cfg.dispInfo(['matRad: Accumulating INTERVAL3 OAR covariance/SVD for ', ...
    '%d voxels in %d batches of up to %d voxels.\n'], ...
    numel(oarRows),numel(batches),batchSize);

for b = 1:numel(batches)
    batchTimer = tic;
    rows = batches{b};
    logBatchProgress(matRad_cfg,cfg,'OAR covariance/SVD',b,numel(batches));
    logBatchStart(matRad_cfg,cfg,'OAR covariance/SVD',b,numel(batches),numel(rows));
    scenarioRows = cell(numel(scenarioDijIx),1);
    for s = 1:numel(scenarioDijIx)
        logScenarioProgress(matRad_cfg,cfg,'OAR covariance/SVD',b,numel(batches), ...
            s,numel(scenarioDijIx),scenarioDijIx(s),scenarioWeights(s));
        scenarioRows{s} = matRad_getDoseIntervalScenarioRows(quantity,scenarioDijIx(s), ...
            scenarioMaps{s},rows,matRad_cfg);
    end

    batchOffset = find(oarRows == rows(1),1,'first') - 1;
    for localIx = 1:numel(rows)
        logVoxelProgress(matRad_cfg,cfg,'OAR covariance/SVD',b,numel(batches), ...
            localIx,numel(rows));
        scenarioMatrixRows = cell(numel(scenarioDijIx),1);
        for s = 1:numel(scenarioDijIx)
            scenarioMatrixRows{s} = scenarioRows{s}(localIx,:);
        end
        scenarioMatrix = vertcat(scenarioMatrixRows{:});

        centerRow = dij_interval.center(rows(localIx),:);
        covariance = scenarioMatrix' * spdiags(scenarioWeights,0, ...
            numel(scenarioWeights),numel(scenarioWeights)) * scenarioMatrix - ...
            centerRow' * centerRow;
        covariance = 0.5 .* (covariance + covariance');

        [U,S,V,k] = truncateCovarianceSvd(covariance,cfg,numBixels);
        intervalIx = batchOffset + localIx;
        dij_interval.k(intervalIx) = k;
        dij_interval.U{intervalIx} = U;
        dij_interval.S{intervalIx} = S;
        dij_interval.V{intervalIx} = V;
    end
    logBatchEnd(matRad_cfg,cfg,'OAR covariance/SVD',b,numel(batches),toc(batchTimer));
end
matRad_cfg.dispInfo('matRad: INTERVAL3 OAR covariance/SVD accumulation finished.\n');
end

function logBatchProgress(matRad_cfg,cfg,stageName,batchIx,numBatches)
if numBatches == 0
    return;
end

progressStep = max(1,ceil(numBatches/10));
if isDetailedProgress(cfg) || batchIx == 1 || batchIx == numBatches || ...
   mod(batchIx,progressStep) == 0
    matRad_cfg.dispInfo('matRad: %s progress %d/%d (%.0f%%).\n', ...
        stageName,batchIx,numBatches,100*batchIx/numBatches);
end
end

function logBatchStart(matRad_cfg,cfg,stageName,batchIx,numBatches,numRows)
if isDetailedProgress(cfg)
    matRad_cfg.dispInfo('matRad: %s batch %d/%d started (%d voxels).\n', ...
        stageName,batchIx,numBatches,numRows);
end
end

function logBatchEnd(matRad_cfg,cfg,stageName,batchIx,numBatches,elapsedTime)
if isDetailedProgress(cfg)
    matRad_cfg.dispInfo('matRad: %s batch %d/%d finished in %.2f s.\n', ...
        stageName,batchIx,numBatches,elapsedTime);
end
end

function logScenarioProgress(matRad_cfg,cfg,stageName,batchIx,numBatches, ...
    scenarioRowIx,numScenarios,scenarioDijIx,scenarioWeight)
if isDetailedProgress(cfg)
    matRad_cfg.dispInfo(['matRad: %s batch %d/%d scenario %d/%d ', ...
        '(linear scenario %d, weight %.6g).\n'],stageName,batchIx, ...
        numBatches,scenarioRowIx,numScenarios,scenarioDijIx,scenarioWeight);
end
end

function logVoxelProgress(matRad_cfg,cfg,stageName,batchIx,numBatches,voxelIx,numVoxels)
if ~isDetailedProgress(cfg) || numVoxels == 0
    return;
end

progressStep = max(1,ceil(numVoxels/10));
if voxelIx == 1 || voxelIx == numVoxels || mod(voxelIx,progressStep) == 0
    matRad_cfg.dispInfo(['matRad: %s batch %d/%d voxel progress ', ...
        '%d/%d (%.0f%%).\n'],stageName,batchIx,numBatches, ...
        voxelIx,numVoxels,100*voxelIx/numVoxels);
end
end

function detailed = isDetailedProgress(cfg)
detailed = isfield(cfg,'ProgressLevel') && strcmp(cfg.ProgressLevel,'detailed');
end

function guardOARSvdMemory(numScenarios,numBixels,cfg,matRad_cfg)
memoryLimitMB = resolveMemoryLimitMB(cfg);
estimatedMB = estimateOARSvdMemoryMB(numScenarios,numBixels,cfg);

matRad_cfg.dispInfo(['matRad: Estimated INTERVAL3 OAR covariance/SVD memory ', ...
    'per voxel is %.2f MB (limit %.2f MB).\n'],estimatedMB,memoryLimitMB);

if estimatedMB > memoryLimitMB
    matRad_cfg.dispError(['INTERVAL3 OAR covariance/SVD estimated memory per ', ...
        'voxel is %.2f MB, which exceeds MemoryLimitMB %.2f MB. Increase ', ...
        'cfg.MemoryLimitMB, reduce cfg.KMax, reduce the number of bixels, ', ...
        'or use INTERVAL2.'],estimatedMB,memoryLimitMB);
end
end

function estimatedMB = estimateOARSvdMemoryMB(numScenarios,numBixels,cfg)
bytesPerDouble = 8;
kMax = min(cfg.KMax,numBixels);

denseCovarianceBytes = 2 * numBixels^2 * bytesPerDouble;
svdWorkspaceBytes = 4 * numBixels^2 * bytesPerDouble;
scenarioMatrixBytes = numScenarios * numBixels * bytesPerDouble;
truncatedSvdBytes = (2*numBixels*kMax + kMax^2) * bytesPerDouble;

estimatedMB = (denseCovarianceBytes + svdWorkspaceBytes + ...
    scenarioMatrixBytes + truncatedSvdBytes) / 1e6;
end

function memoryLimitMB = resolveMemoryLimitMB(cfg)
memoryLimitMB = cfg.MemoryLimitMB;
if isempty(memoryLimitMB)
    memoryLimitMB = 256;
end
end

function batchSize = resolveBatchSize(numRows,numScenarios,numBixels,cfg)
if numRows == 0
    batchSize = 0;
    return;
end

if ~isempty(cfg.BatchSize)
    batchSize = min(numRows,cfg.BatchSize);
    return;
end

memoryLimitMB = resolveMemoryLimitMB(cfg);
bytesPerRow = max(1,numScenarios*numBixels*8);
batchSize = max(1,floor(memoryLimitMB*1e6/bytesPerRow));
batchSize = min(numRows,batchSize);
end

function batches = makeBatches(rows,batchSize)
if isempty(rows)
    batches = {};
    return;
end

numBatches = ceil(numel(rows)/batchSize);
batches = cell(numBatches,1);
for b = 1:numBatches
    firstIx = (b - 1)*batchSize + 1;
    lastIx = min(b*batchSize,numel(rows));
    batches{b} = rows(firstIx:lastIx);
end
end

function [U,S,V,k] = truncateCovarianceSvd(covariance,cfg,numBixels)
tolerance = 1e-12;
if nnz(covariance) == 0 || norm(covariance,'fro') <= tolerance
    [U,S,V,k] = zeroRankSvd(numBixels);
    return;
end

kMax = min(cfg.KMax,numBixels);
if kMax == 0
    [U,S,V,k] = zeroRankSvd(numBixels);
    return;
end

if numBixels <= kMax || numBixels <= 20
    [Ufull,Sfull,Vfull] = svd(full(covariance),'econ');
else
    try
        [Ufull,Sfull,Vfull] = svds(covariance,kMax,'largest');
    catch
        [Ufull,Sfull,Vfull] = svd(full(covariance),'econ');
    end
end

singularValues = diag(Sfull);
positiveIx = find(singularValues > tolerance);
if isempty(positiveIx)
    [U,S,V,k] = zeroRankSvd(numBixels);
    return;
end

Ufull = Ufull(:,positiveIx);
Sfull = Sfull(positiveIx,positiveIx);
Vfull = Vfull(:,positiveIx);
singularValues = singularValues(positiveIx);

if strcmp(cfg.KMode,'dynamic')
    totalEnergy = sum(singularValues.^2);
    k = find(cumsum(singularValues.^2)./totalEnergy >= cfg.RetentionThreshold,1,'first');
    k = min(k,kMax);
else
    k = min(kMax,numel(singularValues));
end

U = sparse(Ufull(:,1:k));
S = sparse(Sfull(1:k,1:k));
V = sparse(Vfull(:,1:k));
end

function [U,S,V,k] = zeroRankSvd(numBixels)
k = 0;
U = sparse(numBixels,0);
S = sparse(0,0);
V = sparse(numBixels,0);
end
