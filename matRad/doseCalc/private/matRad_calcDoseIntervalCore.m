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

if ~isempty(oarRows) && strcmp(intervalMode,'INTERVAL3')
    guardOARSvdMemory(numel(scenarioDijIx),numBixels,cfg,matRad_cfg);
    dij_interval = accumulateOARCenterAndSvd(dij_interval,quantity,scenarioDijIx, ...
        scenarioWeights,scenarioMaps,oarRows,cfg,numBixels,matRad_cfg);
elseif ~isempty(oarRows)
    dij_interval = accumulateOARCenter(dij_interval,quantity,scenarioDijIx, ...
        scenarioWeights,scenarioMaps,oarRows,cfg,numBixels,matRad_cfg);
else
    matRad_cfg.dispInfo('matRad: No OAR voxels selected for interval OAR term.\n');
end

if cfg.CollectTiming
    dij_interval.timing.totalSeconds = toc(timer);
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
if cfg.CollectTiming
    dij_interval.timing = matRad_initializeDoseIntervalTiming('interval', ...
        intervalMode,cfg, ...
        numel(targetRows),numel(oarRows),numel(scenarioDijIx),numBixels);
end

if strcmp(intervalMode,'INTERVAL3')
    dij_interval.k = zeros(numel(oarRows),1);
    dij_interval.U = cell(numel(oarRows),1);
    dij_interval.S = cell(numel(oarRows),1);
    dij_interval.V = cell(numel(oarRows),1);
end
end

function stageTiming = accumulateDoseIntervalBatchTiming(stageTiming,batchTiming)
if isempty(batchTiming)
    return;
end

timingFields = {'extractMapSeconds','centerAccumSeconds', ...
    'radiusMultiplySeconds','svdSeconds','wallSeconds'};
for f = 1:numel(timingFields)
    fieldName = timingFields{f};
    stageTiming.(fieldName) = stageTiming.(fieldName) + ...
        batchTiming.(fieldName);
end
end

function dij_interval = accumulateTargetInterval(dij_interval,quantity,scenarioDijIx, ...
    scenarioWeights,scenarioMaps,targetRows,cfg,numBixels,matRad_cfg)
batchSize = resolveBatchSize(numel(targetRows),numel(scenarioDijIx),numBixels,cfg);
batches = makeBatches(targetRows,batchSize);
numBatches = numel(batches);
matRad_cfg.dispInfo(['matRad: Accumulating target interval center/radius ', ...
    'for %d voxels in %d batches of up to %d voxels.\n'], ...
    numel(targetRows),numBatches,batchSize);
if cfg.CollectTiming
    stageTimer = tic;
    stageTiming = matRad_initializeDoseIntervalTiming('stage','target', ...
        numel(targetRows),batchSize,numBatches);
end

for b = 1:numBatches
    batchTimer = tic;
    rows = batches{b};
    logBatchProgress(matRad_cfg,cfg,'Target interval',b,numBatches);
    logBatchStart(matRad_cfg,cfg,'Target interval',b,numBatches,numel(rows));
    centerBlock = sparse(numel(rows),numBixels);
    radiusBlock = sparse(numBixels,numBixels);

    for s = 1:numel(scenarioDijIx)
        logScenarioProgress(matRad_cfg,cfg,'Target interval',b,numBatches, ...
            s,numel(scenarioDijIx),scenarioDijIx(s),scenarioWeights(s));
        if cfg.CollectTiming
            sectionTimer = tic;
        end
        scenarioRows = matRad_getDoseIntervalScenarioRows(quantity,scenarioDijIx(s), ...
            scenarioMaps{s},rows,matRad_cfg);
        if cfg.CollectTiming
            stageTiming.extractMapSeconds = stageTiming.extractMapSeconds + ...
                toc(sectionTimer);
            sectionTimer = tic;
        end
        centerBlock = centerBlock + scenarioWeights(s) .* scenarioRows;
        if cfg.CollectTiming
            stageTiming.centerAccumSeconds = stageTiming.centerAccumSeconds + ...
                toc(sectionTimer);
            sectionTimer = tic;
        end
        radiusBlock = radiusBlock + scenarioRows' * (scenarioWeights(s) .* scenarioRows);
        if cfg.CollectTiming
            stageTiming.radiusMultiplySeconds = ...
                stageTiming.radiusMultiplySeconds + toc(sectionTimer);
        end
    end

    dij_interval.center(rows,:) = centerBlock;
    dij_interval.radius = dij_interval.radius + radiusBlock;
    logBatchEnd(matRad_cfg,cfg,'Target interval',b,numBatches,toc(batchTimer));
end
if cfg.CollectTiming
    stageTiming.wallSeconds = toc(stageTimer);
    dij_interval.timing.target = stageTiming;
end
matRad_cfg.dispInfo('matRad: Target interval accumulation finished.\n');
end

function dij_interval = accumulateOARCenter(dij_interval,quantity,scenarioDijIx, ...
    scenarioWeights,scenarioMaps,oarRows,cfg,numBixels,matRad_cfg)
batchSize = resolveBatchSize(numel(oarRows),numel(scenarioDijIx),numBixels,cfg);
batches = makeBatches(oarRows,batchSize);
numBatches = numel(batches);
matRad_cfg.dispInfo(['matRad: Accumulating OAR interval center for %d voxels ', ...
    'in %d batches of up to %d voxels.\n'], ...
    numel(oarRows),numBatches,batchSize);
if cfg.CollectTiming
    stageTimer = tic;
    stageTiming = matRad_initializeDoseIntervalTiming('stage','oar', ...
        numel(oarRows),batchSize,numBatches);
end

for b = 1:numBatches
    batchTimer = tic;
    rows = batches{b};
    logBatchProgress(matRad_cfg,cfg,'OAR center',b,numBatches);
    logBatchStart(matRad_cfg,cfg,'OAR center',b,numBatches,numel(rows));
    centerBlock = sparse(numel(rows),numBixels);
    for s = 1:numel(scenarioDijIx)
        logScenarioProgress(matRad_cfg,cfg,'OAR center',b,numBatches, ...
            s,numel(scenarioDijIx),scenarioDijIx(s),scenarioWeights(s));
        if cfg.CollectTiming
            sectionTimer = tic;
        end
        scenarioRows = matRad_getDoseIntervalScenarioRows(quantity,scenarioDijIx(s), ...
            scenarioMaps{s},rows,matRad_cfg);
        if cfg.CollectTiming
            stageTiming.extractMapSeconds = stageTiming.extractMapSeconds + ...
                toc(sectionTimer);
            sectionTimer = tic;
        end
        centerBlock = centerBlock + scenarioWeights(s) .* scenarioRows;
        if cfg.CollectTiming
            stageTiming.centerAccumSeconds = stageTiming.centerAccumSeconds + ...
                toc(sectionTimer);
        end
    end
    dij_interval.center(rows,:) = centerBlock;
    logBatchEnd(matRad_cfg,cfg,'OAR center',b,numBatches,toc(batchTimer));
end
if cfg.CollectTiming
    stageTiming.wallSeconds = toc(stageTimer);
    dij_interval.timing.oar = stageTiming;
end
matRad_cfg.dispInfo('matRad: OAR interval center accumulation finished.\n');
end

function dij_interval = accumulateOARCenterAndSvd(dij_interval,quantity,scenarioDijIx, ...
    scenarioWeights,scenarioMaps,oarRows,cfg,numBixels,matRad_cfg)
batchSize = resolveBatchSize(numel(oarRows),numel(scenarioDijIx),numBixels,cfg);
batches = makeBatches(oarRows,batchSize);
numBatches = numel(batches);
matRad_cfg.dispInfo(['matRad: Accumulating INTERVAL3 OAR interval center ', ...
    'and scenario SVD for %d voxels in %d batches of up to %d voxels.\n'], ...
    numel(oarRows),numBatches,batchSize);

if cfg.CollectTiming
    stageTimer = tic;
    stageTiming = matRad_initializeDoseIntervalTiming('stage','oar', ...
        numel(oarRows),batchSize,numBatches);
    sectionTimer = tic;
end
useParallel = configureOARSvdParallel(numBatches,batchSize,quantity, ...
    scenarioMaps,numel(scenarioDijIx),numBixels,cfg,matRad_cfg);
if cfg.CollectTiming
    stageTiming.parallelSetupSeconds = toc(sectionTimer);
    stageTiming.parallelEnabled = useParallel;
end

if useParallel
    centerBlocks = cell(numBatches,1);
    kBlocks = cell(numBatches,1);
    UBlocks = cell(numBatches,1);
    SBlocks = cell(numBatches,1);
    VBlocks = cell(numBatches,1);
    timingBlocks = cell(numBatches,1);

    logLevel = matRad_cfg.logLevel;
    progressQueue = [];
    progressListener = [];
    nFinishedBatches = 0;
    if exist('parallel.pool.DataQueue','class') == 8
        progressQueue = parallel.pool.DataQueue;
        progressListener = afterEach(progressQueue,@(~) updateParallelProgress());
    end

    if cfg.CollectTiming
        sectionTimer = tic;
    end
    parfor b = 1:numBatches
        workerCfg = MatRad_Config.instance();
        workerCfg.logLevel = logLevel;

        [centerBlock,kBlock,UBlock,SBlock,VBlock,batchTiming] = ...
            computeOARCenterAndSvdBatch( ...
            quantity,scenarioDijIx,scenarioWeights,scenarioMaps,batches{b}, ...
            cfg,numBixels,workerCfg,'OAR center/scenario SVD',b, ...
            numBatches,false);
        centerBlocks{b} = centerBlock;
        kBlocks{b} = kBlock;
        UBlocks{b} = UBlock;
        SBlocks{b} = SBlock;
        VBlocks{b} = VBlock;
        timingBlocks{b} = batchTiming;

        if ~isempty(progressQueue)
            send(progressQueue,b);
        end
    end
    if cfg.CollectTiming
        stageTiming.parallelComputeWallSeconds = toc(sectionTimer);
    end

    if ~isempty(progressListener)
        delete(progressListener);
    end

    if cfg.CollectTiming
        sectionTimer = tic;
    end
    dij_interval = assembleOARCenterAndSvdBatches(dij_interval,oarRows,batches, ...
        centerBlocks,kBlocks,UBlocks,SBlocks,VBlocks);
    if cfg.CollectTiming
        stageTiming.serialAssemblySeconds = toc(sectionTimer);
        for b = 1:numBatches
            stageTiming = accumulateDoseIntervalBatchTiming(stageTiming, ...
                timingBlocks{b});
        end
    end
else
    centerBlocks = cell(numBatches,1);
    kBlocks = cell(numBatches,1);
    UBlocks = cell(numBatches,1);
    SBlocks = cell(numBatches,1);
    VBlocks = cell(numBatches,1);
    timingBlocks = cell(numBatches,1);

    for b = 1:numBatches
        batchTimer = tic;
        rows = batches{b};
        logBatchProgress(matRad_cfg,cfg,'OAR center/scenario SVD',b,numBatches);
        logBatchStart(matRad_cfg,cfg,'OAR center/scenario SVD',b,numBatches,numel(rows));
        [centerBlock,kBlock,UBlock,SBlock,VBlock,batchTiming] = ...
            computeOARCenterAndSvdBatch( ...
            quantity,scenarioDijIx,scenarioWeights,scenarioMaps,rows,cfg, ...
            numBixels,matRad_cfg,'OAR center/scenario SVD',b,numBatches,true);
        centerBlocks{b} = centerBlock;
        kBlocks{b} = kBlock;
        UBlocks{b} = UBlock;
        SBlocks{b} = SBlock;
        VBlocks{b} = VBlock;
        timingBlocks{b} = batchTiming;
        logBatchEnd(matRad_cfg,cfg,'OAR center/scenario SVD',b,numBatches, ...
            toc(batchTimer));
    end

    if cfg.CollectTiming
        sectionTimer = tic;
    end
    dij_interval = assembleOARCenterAndSvdBatches(dij_interval,oarRows,batches, ...
        centerBlocks,kBlocks,UBlocks,SBlocks,VBlocks);
    if cfg.CollectTiming
        stageTiming.serialAssemblySeconds = toc(sectionTimer);
        for b = 1:numBatches
            stageTiming = accumulateDoseIntervalBatchTiming(stageTiming, ...
                timingBlocks{b});
        end
    end
end

if cfg.CollectTiming
    stageTiming.wallSeconds = toc(stageTimer);
    dij_interval.timing.oar = stageTiming;
end

matRad_cfg.dispInfo(['matRad: INTERVAL3 OAR interval center and scenario ', ...
    'SVD accumulation finished.\n']);

    function updateParallelProgress()
        nFinishedBatches = nFinishedBatches + 1;
        logBatchProgress(matRad_cfg,cfg,'OAR center/scenario SVD', ...
            nFinishedBatches,numBatches);
        drawnow('limitrate');
    end
end

function [centerBlock,kBlock,UBlock,SBlock,VBlock,batchTiming] = ...
    computeOARCenterAndSvdBatch( ...
    quantity,scenarioDijIx,scenarioWeights,scenarioMaps,rows,cfg,numBixels, ...
    matRad_cfg,stageName,batchIx,numBatches,enableProgress)
if cfg.CollectTiming
    batchTimer = tic;
    batchTiming = matRad_initializeDoseIntervalTiming('batch',numel(rows));
else
    batchTiming = [];
end
centerBlock = sparse(numel(rows),numBixels);
scenarioRows = cell(numel(scenarioDijIx),1);
for s = 1:numel(scenarioDijIx)
    if enableProgress
        logScenarioProgress(matRad_cfg,cfg,stageName,batchIx,numBatches, ...
            s,numel(scenarioDijIx),scenarioDijIx(s),scenarioWeights(s));
    end
    if cfg.CollectTiming
        sectionTimer = tic;
    end
    scenarioRows{s} = matRad_getDoseIntervalScenarioRows(quantity,scenarioDijIx(s), ...
        scenarioMaps{s},rows,matRad_cfg);
    if cfg.CollectTiming
        batchTiming.extractMapSeconds = batchTiming.extractMapSeconds + ...
            toc(sectionTimer);
        sectionTimer = tic;
    end
    centerBlock = centerBlock + scenarioWeights(s) .* scenarioRows{s};
    if cfg.CollectTiming
        batchTiming.centerAccumSeconds = batchTiming.centerAccumSeconds + ...
            toc(sectionTimer);
    end
end

kBlock = zeros(numel(rows),1);
UBlock = cell(numel(rows),1);
SBlock = cell(numel(rows),1);
VBlock = cell(numel(rows),1);
for localIx = 1:numel(rows)
    if enableProgress
        logVoxelProgress(matRad_cfg,cfg,stageName,batchIx,numBatches, ...
            localIx,numel(rows));
    end

    scenarioMatrixRows = cell(numel(scenarioDijIx),1);
    for s = 1:numel(scenarioDijIx)
        scenarioMatrixRows{s} = scenarioRows{s}(localIx,:);
    end
    scenarioMatrix = vertcat(scenarioMatrixRows{:});

    centerRow = centerBlock(localIx,:);
    if cfg.CollectTiming
        sectionTimer = tic;
    end
    [U,S,V,k] = truncateWeightedScenarioSvd(scenarioMatrix,centerRow, ...
        scenarioWeights,cfg,numBixels);
    if cfg.CollectTiming
        batchTiming.svdSeconds = batchTiming.svdSeconds + toc(sectionTimer);
    end
    kBlock(localIx) = k;
    UBlock{localIx} = U;
    SBlock{localIx} = S;
    VBlock{localIx} = V;
end
if cfg.CollectTiming
    batchTiming.wallSeconds = toc(batchTimer);
end
end

function dij_interval = assembleOARCenterAndSvdBatches(dij_interval,oarRows, ...
    batches,centerBlocks,kBlocks,UBlocks,SBlocks,VBlocks)
intervalOffset = 0;
for b = 1:numel(batches)
    rows = batches{b};
    numRows = numel(rows);
    intervalIx = intervalOffset + (1:numRows);
    dij_interval.center(rows,:) = centerBlocks{b};
    dij_interval.k(intervalIx) = kBlocks{b};
    dij_interval.U(intervalIx) = UBlocks{b};
    dij_interval.S(intervalIx) = SBlocks{b};
    dij_interval.V(intervalIx) = VBlocks{b};
    intervalOffset = intervalOffset + numRows;
end

if intervalOffset ~= numel(oarRows)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('OAR SVD batch assembly did not cover all OAR voxels.');
end
end

function useParallel = configureOARSvdParallel(numBatches,batchSize,quantity, ...
    scenarioMaps,numScenarios,numBixels,cfg,matRad_cfg)
useParallel = false;
if ~cfg.UseParallel || numBatches < 2
    return;
end

if ~isParallelComputingAvailable()
    matRad_cfg.dispWarning(['UseParallel was requested for INTERVAL3 OAR ', ...
        'scenario SVD, but the Parallel Computing Toolbox is unavailable. ', ...
        'Falling back to serial batch accumulation.\n']);
    return;
end

workerMemoryBytes = estimateOARSvdBatchWorkerMemoryBytes(numScenarios, ...
    numBixels,batchSize,cfg) + estimateOARSvdParallelBroadcastBytes( ...
    quantity,scenarioMaps);
[poolSizeLimit,memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount( ...
    workerMemoryBytes,'numTasks',numBatches,'safetyFactor',1.2, ...
    'minWorkerMemoryBytes',512 * 1024^2);

matRad_cfg.dispInfo(['matRad: Estimated INTERVAL3 OAR scenario SVD memory ', ...
    'per parallel worker is %.2f MB, including broadcast data.\n'], ...
    memoryEstimate.workerBytes / 1e6);

if isempty(poolSizeLimit) || poolSizeLimit < 2
    matRad_cfg.dispWarning(['UseParallel was requested for INTERVAL3 OAR ', ...
        'scenario SVD, but the estimated memory only allows one worker. ', ...
        'Falling back to serial batch accumulation.\n']);
    return;
end

try
    pPool = configureDoseIntervalParallelPool(poolSizeLimit,matRad_cfg);
catch
    matRad_cfg.dispWarning(['Could not configure a parallel pool for ', ...
        'INTERVAL3 OAR scenario SVD. Falling back to serial batch ', ...
        'accumulation.\n']);
    return;
end

if isempty(pPool) || pPool.NumWorkers < 2
    matRad_cfg.dispWarning(['UseParallel was requested for INTERVAL3 OAR ', ...
        'scenario SVD, but fewer than two workers are available. Falling ', ...
        'back to serial batch accumulation.\n']);
    return;
end

matRad_cfg.dispInfo(['matRad: Parallel INTERVAL3 OAR scenario SVD uses %d ', ...
    'worker(s) for %d batches.\n'],pPool.NumWorkers,numBatches);
useParallel = true;
end

function available = isParallelComputingAvailable()
available = false;
if exist('parpool','file') ~= 2 || exist('gcp','file') ~= 2
    return;
end

try
    [available,~] = license('checkout','Distrib_Computing_Toolbox');
catch
    available = false;
end

if isempty(available)
    available = false;
end
available = logical(available);
end

function pPool = configureDoseIntervalParallelPool(poolSizeLimit,matRad_cfg)
pPool = gcp('nocreate');

if isempty(pPool)
    pPool = parpool(poolSizeLimit);
elseif pPool.NumWorkers > poolSizeLimit
    matRad_cfg.dispWarning(['Reducing parallel pool from ', ...
        num2str(pPool.NumWorkers),' to ',num2str(poolSizeLimit), ...
        ' worker(s) to keep INTERVAL3 OAR scenario SVD within the ', ...
        'estimated available memory.\n']);
    delete(pPool);
    pPool = parpool(poolSizeLimit);
end
end

function workerMemoryBytes = estimateOARSvdBatchWorkerMemoryBytes(numScenarios, ...
    numBixels,batchSize,cfg)
bytesPerDouble = 8;
kMax = min([cfg.KMax,numBixels,numScenarios]);

scenarioRowsBytes = numScenarios * batchSize * numBixels * bytesPerDouble;
centerBlockBytes = batchSize * numBixels * bytesPerDouble;
perVoxelWorkspaceBytes = estimateOARSvdMemoryMB(numScenarios,numBixels,cfg) * 1e6;
batchFactorBytes = batchSize * (3*numBixels*kMax + kMax^2) * bytesPerDouble;

workerMemoryBytes = scenarioRowsBytes + centerBlockBytes + ...
    perVoxelWorkspaceBytes + batchFactorBytes;
end

function broadcastBytes = estimateOARSvdParallelBroadcastBytes(quantity,scenarioMaps)
broadcastBytes = estimateVariableBytes(quantity) + estimateVariableBytes(scenarioMaps);
end

function bytes = estimateVariableBytes(value)
variableInfo = whos('value');
bytes = variableInfo.bytes;
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

matRad_cfg.dispInfo(['matRad: Estimated INTERVAL3 OAR scenario SVD memory ', ...
    'per voxel is %.2f MB (limit %.2f MB).\n'],estimatedMB,memoryLimitMB);

if estimatedMB > memoryLimitMB
    matRad_cfg.dispError(['INTERVAL3 OAR scenario SVD estimated memory per ', ...
        'voxel is %.2f MB, which exceeds MemoryLimitMB %.2f MB. Increase ', ...
        'cfg.MemoryLimitMB, reduce cfg.KMax, reduce the number of scenarios ', ...
        'or bixels, or use INTERVAL2.'],estimatedMB,memoryLimitMB);
end
end

function estimatedMB = estimateOARSvdMemoryMB(numScenarios,numBixels,cfg)
bytesPerDouble = 8;
kMax = min([cfg.KMax,numBixels,numScenarios]);

scenarioMatrixBytes = 3 * numScenarios * numBixels * bytesPerDouble;
scenarioGramBytes = 3 * numScenarios^2 * bytesPerDouble;
truncatedSvdBytes = (3*numBixels*kMax + kMax^2) * bytesPerDouble;

estimatedMB = (scenarioMatrixBytes + scenarioGramBytes + ...
    truncatedSvdBytes) / 1e6;
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

function [U,S,V,k] = truncateWeightedScenarioSvd(scenarioMatrix,centerRow, ...
    scenarioWeights,cfg,numBixels)
tolerance = 1e-12;

kMax = min(cfg.KMax,numBixels);
if kMax == 0
    [U,S,V,k] = zeroRankSvd(numBixels);
    return;
end

positiveWeightIx = scenarioWeights(:) > 0;
scenarioWeights = scenarioWeights(positiveWeightIx);
scenarioMatrix = scenarioMatrix(positiveWeightIx,:);

activeBixels = unique([find(any(scenarioMatrix,1)) find(centerRow)]);
if isempty(activeBixels)
    [U,S,V,k] = zeroRankSvd(numBixels);
    return;
end

scenarioMatrix = scenarioMatrix(:,activeBixels);
centerRow = centerRow(:,activeBixels);
numScenarios = size(scenarioMatrix,1);

centeredScenarioMatrix = scenarioMatrix - repmat(centerRow,numScenarios,1);
weightedScenarioMatrix = spdiags(sqrt(scenarioWeights),0, ...
    numScenarios,numScenarios) * centeredScenarioMatrix;

if nnz(weightedScenarioMatrix) == 0 || ...
        norm(weightedScenarioMatrix,'fro') <= tolerance
    [U,S,V,k] = zeroRankSvd(numBixels);
    return;
end

% X'*X and X*X' share non-zero eigenvalues; diagonalize the smaller
% scenario Gram and expand right singular vectors back to bixel space.
scenarioGram = weightedScenarioMatrix * weightedScenarioMatrix';
scenarioGram = full(0.5 .* (scenarioGram + scenarioGram'));
[leftSingularVectors,eigenValueMatrix] = eig(scenarioGram);
singularValues = real(diag(eigenValueMatrix));
singularValues(abs(singularValues) <= tolerance) = 0;
singularValues(singularValues < 0) = 0;
[singularValues,sortIx] = sort(singularValues,'descend');
leftSingularVectors = leftSingularVectors(:,sortIx);

positiveIx = find(singularValues > tolerance);
if isempty(positiveIx)
    [U,S,V,k] = zeroRankSvd(numBixels);
    return;
end

singularValues = singularValues(positiveIx);
leftSingularVectors = leftSingularVectors(:,positiveIx);

if strcmp(cfg.KMode,'dynamic')
    k = selectScenarioSvdRank(singularValues,kMax, ...
        cfg.RetentionThreshold,tolerance);
else
    k = min(kMax,numel(singularValues));
end

if k == 0
    [U,S,V,k] = zeroRankSvd(numBixels);
    return;
end

singularValues = singularValues(1:k);
leftSingularVectors = leftSingularVectors(:,1:k);
rightSingularVectors = weightedScenarioMatrix' * leftSingularVectors;
rightSingularVectors = bsxfun(@rdivide,rightSingularVectors, ...
    sqrt(singularValues(:))');

V = sparse(double(numBixels),double(k));
V(activeBixels,:) = sparse(rightSingularVectors);
U = V;
S = spdiags(singularValues(:),0,k,k);
end

function k = selectScenarioSvdRank(singularValues,kMax,retentionThreshold,tolerance)
maxRank = min(kMax,numel(singularValues));
if maxRank == 0
    k = 0;
    return;
end

rankSingularValues = singularValues(1:maxRank);
energyScale = max(abs(rankSingularValues));
if ~isfinite(energyScale) || energyScale <= tolerance
    k = maxRank;
    return;
end

relativeEnergy = (rankSingularValues ./ energyScale).^2;
totalEnergy = sum(relativeEnergy);
if ~isfinite(totalEnergy) || totalEnergy <= 0
    k = maxRank;
    return;
end

thresholdEnergy = retentionThreshold * totalEnergy;
k = find(cumsum(relativeEnergy) >= thresholdEnergy,1,'first');
if isempty(k)
    k = maxRank;
end
end

function [U,S,V,k] = zeroRankSvd(numBixels)
k = 0;
U = sparse(numBixels,0);
S = sparse(0,0);
V = sparse(numBixels,0);
end
