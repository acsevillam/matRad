function [pln_interval,dij_intervalContext] = matRad_calcDoseInterval(ct,cst,stf,pln,dij,cfg)
% matRad_calcDoseInterval calculates generic in-memory interval dose influence data
%
% call
%   [pln_interval,dij_intervalContext] = matRad_calcDoseInterval(ct,cst,stf,pln,dij,cfg)
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
%   cfg:          dose interval configuration struct with required field
%                 IntervalMode: 'INTERVAL2' or 'INTERVAL3'
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

if nargin < 6 || isempty(cfg)
    cfg = struct();
end
matRad_cfg = MatRad_Config.instance();
if ~isstruct(cfg)
    matRad_cfg.dispError('Dose interval configuration must be a struct.');
end
if ~isfield(cfg,'IntervalMode') || isempty(cfg.IntervalMode)
    matRad_cfg.dispError('cfg.IntervalMode must be ''INTERVAL2'' or ''INTERVAL3''.');
end
intervalMode = cfg.IntervalMode;
intervalMode = normalizeIntervalMode(intervalMode,matRad_cfg);
cfg.IntervalMode = intervalMode;

timer = tic;

ctx = matRad_resolveScenarioDoseInputs(ct,cst,pln,dij,cfg, ...
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
matRad_cfg.dispInfo('matRad: Dose interval radius mode is ''%s''.\n', ...
    cfg.RadiusMode);
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
    guardOARRadiusFactorMemory(numel(scenarioDijIx),numBixels,cfg,matRad_cfg);
    dij_interval = accumulateOARCenterAndRadiusFactor(dij_interval,quantity,scenarioDijIx, ...
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
dij_intervalContext = buildIntervalDijContext(ct,cst,stf,pln,dij, ...
    dij_interval,quantity,cfg);
pln_interval.multScen = dij_intervalContext.scenarioModel;

matRad_cfg.dispInfo('matRad: Finished %s dose interval calculation in %.2f s.\n', ...
    intervalMode,toc(timer));

end

function intervalMode = normalizeIntervalMode(intervalMode,matRad_cfg)
if isstring(intervalMode) && isscalar(intervalMode)
    intervalMode = char(intervalMode);
end
if ~ischar(intervalMode) || isempty(intervalMode)
    matRad_cfg.dispError('intervalMode must be ''INTERVAL2'' or ''INTERVAL3''.');
end

intervalMode = upper(intervalMode);
if ~any(strcmp(intervalMode,{'INTERVAL2','INTERVAL3'}))
    matRad_cfg.dispError('intervalMode must be ''INTERVAL2'' or ''INTERVAL3''.');
end
end

function dij_intervalContext = buildIntervalDijContext(ct,cst,stf,pln,dij, ...
    dij_interval,quantity,cfg)
numBixels = size(dij_interval.center,2);
dij_intervalContext = matRad_buildNominalDijContext( ...
    ct,cst,stf,pln,cfg,numBixels,dij);
dij_intervalContext.intervalQuantity = quantity.name;
dij_intervalContext.intervalQuantityField = quantity.field;
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
dij_interval.radiusMode = cfg.RadiusMode;
if cfg.CollectTiming
    dij_interval.timing = matRad_initializeDoseIntervalTiming('interval', ...
        intervalMode,cfg, ...
        numel(targetRows),numel(oarRows),numel(scenarioDijIx),numBixels);
end

if strcmp(intervalMode,'INTERVAL3')
    dij_interval.OARRadiusRank = zeros(numel(oarRows),1);
    dij_interval.OARRadiusFactor = cell(numel(oarRows),1);
end
end

function stageTiming = accumulateDoseIntervalBatchTiming(stageTiming,batchTiming)
if isempty(batchTiming)
    return;
end

timingFields = {'extractMapSeconds','centerAccumSeconds', ...
    'radiusMultiplySeconds','factorSeconds','centeredRowsSeconds', ...
    'wallSeconds'};
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
    [~,meanRows,batchStatsTiming,secondMoment,extremeDeltaRows] = computeIntervalScenarioDoseBatchStats( ...
        quantity,scenarioDijIx,scenarioWeights,scenarioMaps,rows,numBixels, ...
        cfg,matRad_cfg,'Target interval',b,numBatches,true,true);
    if cfg.CollectTiming
        stageTiming.extractMapSeconds = stageTiming.extractMapSeconds + ...
            batchStatsTiming.extractMapSeconds;
        stageTiming.centerAccumSeconds = stageTiming.centerAccumSeconds + ...
            batchStatsTiming.centerAccumSeconds;
        stageTiming.radiusMultiplySeconds = stageTiming.radiusMultiplySeconds + ...
            batchStatsTiming.radiusMultiplySeconds;
        sectionTimer = tic;
    end

    radiusBlock = computeTargetRadiusBlock(meanRows,secondMoment, ...
        extremeDeltaRows,cfg,numBixels);
    if cfg.CollectTiming
        stageTiming.radiusMultiplySeconds = ...
            stageTiming.radiusMultiplySeconds + toc(sectionTimer);
    end

    dij_interval.center(rows,:) = meanRows;
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
    [~,meanRows,batchStatsTiming] = computeIntervalScenarioDoseBatchStats( ...
        quantity,scenarioDijIx,scenarioWeights,scenarioMaps,rows,numBixels, ...
        cfg,matRad_cfg,'OAR center',b,numBatches,true,false);
    if cfg.CollectTiming
        stageTiming.extractMapSeconds = stageTiming.extractMapSeconds + ...
            batchStatsTiming.extractMapSeconds;
        stageTiming.centerAccumSeconds = stageTiming.centerAccumSeconds + ...
            batchStatsTiming.centerAccumSeconds;
    end
    dij_interval.center(rows,:) = meanRows;
    logBatchEnd(matRad_cfg,cfg,'OAR center',b,numBatches,toc(batchTimer));
end
if cfg.CollectTiming
    stageTiming.wallSeconds = toc(stageTimer);
    dij_interval.timing.oar = stageTiming;
end
matRad_cfg.dispInfo('matRad: OAR interval center accumulation finished.\n');
end

function [scenarioRows,meanRows,batchTiming,secondMoment,extremeDeltaRows,centeredRows] = computeIntervalScenarioDoseBatchStats( ...
    quantity,scenarioDijIx,scenarioWeights,scenarioMaps,rows,numBixels, ...
    cfg,matRad_cfg,stageName,batchIx,numBatches,enableProgress, ...
    computeRadiusStats,computeCenteredRows)
if nargin < 14
    computeCenteredRows = false;
end

if cfg.CollectTiming
    batchTiming = matRad_initializeDoseIntervalTiming('batch',numel(rows));
else
    batchTiming = [];
end

options = struct();
options.StoreScenarioRows = true;
options.ComputeSecondMoment = computeRadiusStats && strcmp(cfg.RadiusMode,'std');
options.ComputeCenteredRows = computeCenteredRows;
options.ComputeExtremeDelta = computeRadiusStats && strcmp(cfg.RadiusMode,'extreme');
options.CollectTiming = cfg.CollectTiming;
if enableProgress
    options.ScenarioProgressFcn = @(s) logScenarioProgress(matRad_cfg,cfg, ...
        stageName,batchIx,numBatches,s,numel(scenarioDijIx), ...
        scenarioDijIx(s),scenarioWeights(s));
end

[stats,helperTiming] = matRad_computeScenarioDoseBatchStats(quantity, ...
    scenarioDijIx,scenarioWeights,scenarioMaps,rows,numBixels, ...
    matRad_cfg,options);

scenarioRows = stats.scenarioRows;
meanRows = stats.meanRows;
secondMoment = stats.secondMoment;
extremeDeltaRows = stats.extremeDeltaRows;
centeredRows = stats.centeredRows;

if cfg.CollectTiming
    batchTiming.extractMapSeconds = helperTiming.extractMapSeconds;
    batchTiming.centerAccumSeconds = helperTiming.centerAccumSeconds;
    batchTiming.radiusMultiplySeconds = helperTiming.secondMomentSeconds + ...
        helperTiming.extremeDeltaSeconds;
    batchTiming.centeredRowsSeconds = helperTiming.centeredRowsSeconds;
end
end

function radiusBlock = computeTargetRadiusBlock(meanRows,secondMoment, ...
    extremeDeltaRows,cfg,numBixels)
if strcmp(cfg.RadiusMode,'std')
    radiusBlock = secondMoment;
    return;
end

deltaSquared = full(sum(extremeDeltaRows.^2,1));
radiusBlock = meanRows' * meanRows + ...
    spdiags(deltaSquared(:),0,numBixels,numBixels);
end

function dij_interval = accumulateOARCenterAndRadiusFactor(dij_interval,quantity,scenarioDijIx, ...
    scenarioWeights,scenarioMaps,oarRows,cfg,numBixels,matRad_cfg)
batchSize = resolveBatchSize(numel(oarRows),numel(scenarioDijIx),numBixels,cfg);
batches = makeBatches(oarRows,batchSize);
numBatches = numel(batches);
matRad_cfg.dispInfo(['matRad: Accumulating INTERVAL3 OAR interval center ', ...
    'and radius factor for %d voxels in %d batches of up to %d voxels.\n'], ...
    numel(oarRows),numBatches,batchSize);

if cfg.CollectTiming
    stageTimer = tic;
    stageTiming = matRad_initializeDoseIntervalTiming('stage','oar', ...
        numel(oarRows),batchSize,numBatches);
    sectionTimer = tic;
end
useParallel = configureOARRadiusFactorParallel(numBatches,batchSize,quantity, ...
    scenarioMaps,numel(scenarioDijIx),numBixels,cfg,matRad_cfg);
if cfg.CollectTiming
    stageTiming.parallelSetupSeconds = toc(sectionTimer);
    stageTiming.parallelEnabled = useParallel;
end

if useParallel
    centerBlocks = cell(numBatches,1);
    rankBlocks = cell(numBatches,1);
    factorBlocks = cell(numBatches,1);
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

        [centerBlock,rankBlock,factorBlock,batchTiming] = ...
            computeOARCenterAndRadiusFactorBatch( ...
            quantity,scenarioDijIx,scenarioWeights,scenarioMaps,batches{b}, ...
            cfg,numBixels,workerCfg,'OAR center/radius factor',b, ...
            numBatches,false);
        centerBlocks{b} = centerBlock;
        rankBlocks{b} = rankBlock;
        factorBlocks{b} = factorBlock;
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
    dij_interval = assembleOARCenterAndRadiusFactorBatches(dij_interval,oarRows,batches, ...
        centerBlocks,rankBlocks,factorBlocks);
    if cfg.CollectTiming
        stageTiming.serialAssemblySeconds = toc(sectionTimer);
        for b = 1:numBatches
            stageTiming = accumulateDoseIntervalBatchTiming(stageTiming, ...
                timingBlocks{b});
        end
    end
else
    centerBlocks = cell(numBatches,1);
    rankBlocks = cell(numBatches,1);
    factorBlocks = cell(numBatches,1);
    timingBlocks = cell(numBatches,1);

    for b = 1:numBatches
        batchTimer = tic;
        rows = batches{b};
        logBatchProgress(matRad_cfg,cfg,'OAR center/radius factor',b,numBatches);
        logBatchStart(matRad_cfg,cfg,'OAR center/radius factor',b,numBatches,numel(rows));
        [centerBlock,rankBlock,factorBlock,batchTiming] = ...
            computeOARCenterAndRadiusFactorBatch( ...
            quantity,scenarioDijIx,scenarioWeights,scenarioMaps,rows,cfg, ...
            numBixels,matRad_cfg,'OAR center/radius factor',b,numBatches,true);
        centerBlocks{b} = centerBlock;
        rankBlocks{b} = rankBlock;
        factorBlocks{b} = factorBlock;
        timingBlocks{b} = batchTiming;
        logBatchEnd(matRad_cfg,cfg,'OAR center/radius factor',b,numBatches, ...
            toc(batchTimer));
    end

    if cfg.CollectTiming
        sectionTimer = tic;
    end
    dij_interval = assembleOARCenterAndRadiusFactorBatches(dij_interval,oarRows,batches, ...
        centerBlocks,rankBlocks,factorBlocks);
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
    'radius factor accumulation finished.\n']);

    function updateParallelProgress()
        nFinishedBatches = nFinishedBatches + 1;
        logBatchProgress(matRad_cfg,cfg,'OAR center/radius factor', ...
            nFinishedBatches,numBatches);
        drawnow('limitrate');
    end
end

function [centerBlock,rankBlock,factorBlock,batchTiming] = ...
    computeOARCenterAndRadiusFactorBatch( ...
    quantity,scenarioDijIx,scenarioWeights,scenarioMaps,rows,cfg,numBixels, ...
    matRad_cfg,stageName,batchIx,numBatches,enableProgress)
if cfg.CollectTiming
    batchTimer = tic;
end
[~,centerBlock,batchTiming,~,~,centeredRows] = computeIntervalScenarioDoseBatchStats( ...
    quantity,scenarioDijIx,scenarioWeights,scenarioMaps,rows,numBixels, ...
    cfg,matRad_cfg,stageName,batchIx,numBatches,enableProgress,false,true);

rankBlock = zeros(numel(rows),1);
factorBlock = cell(numel(rows),1);
for localIx = 1:numel(rows)
    if enableProgress
        logVoxelProgress(matRad_cfg,cfg,stageName,batchIx,numBatches, ...
            localIx,numel(rows));
    end

    centeredMatrixRows = cell(numel(scenarioDijIx),1);
    for s = 1:numel(scenarioDijIx)
        centeredMatrixRows{s} = centeredRows{s}(localIx,:);
    end
    centeredScenarioMatrix = vertcat(centeredMatrixRows{:});

    if cfg.CollectTiming
        sectionTimer = tic;
    end
    [factor,rank] = buildWeightedScenarioRadiusFactor( ...
        centeredScenarioMatrix,scenarioWeights,cfg,numBixels);
    if cfg.CollectTiming
        batchTiming.factorSeconds = batchTiming.factorSeconds + ...
            toc(sectionTimer);
    end
    rankBlock(localIx) = rank;
    factorBlock{localIx} = factor;
end

if cfg.CollectTiming
    batchTiming.wallSeconds = toc(batchTimer);
end
end

function dij_interval = assembleOARCenterAndRadiusFactorBatches(dij_interval,oarRows, ...
    batches,centerBlocks,rankBlocks,factorBlocks)
intervalOffset = 0;
for b = 1:numel(batches)
    rows = batches{b};
    numRows = numel(rows);
    intervalIx = intervalOffset + (1:numRows);
    dij_interval.center(rows,:) = centerBlocks{b};
    dij_interval.OARRadiusRank(intervalIx) = rankBlocks{b};
    dij_interval.OARRadiusFactor(intervalIx) = factorBlocks{b};
    intervalOffset = intervalOffset + numRows;
end

if intervalOffset ~= numel(oarRows)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('OAR radius factor batch assembly did not cover all OAR voxels.');
end
end

function useParallel = configureOARRadiusFactorParallel(numBatches,batchSize,quantity, ...
    scenarioMaps,numScenarios,numBixels,cfg,matRad_cfg)
useParallel = false;
if ~cfg.UseParallel || numBatches < 2
    return;
end

workerMemoryBytes = estimateOARRadiusFactorBatchWorkerMemoryBytes(numScenarios, ...
    numBixels,batchSize,cfg) + estimateOARRadiusFactorParallelBroadcastBytes( ...
    quantity,scenarioMaps);
parallelOptions = matRad_doseParallelPoolOptions( ...
    cfg,matRad_cfg,'parallelOptions');
[useParallel,~,~] = matRad_configureSafeDoseParallelPool( ...
    workerMemoryBytes,numBatches,matRad_cfg, ...
    'INTERVAL3 OAR radius factor', ...
    'fallbackDescription','serial batch accumulation', ...
    parallelOptions{:});
end

function workerMemoryBytes = estimateOARRadiusFactorBatchWorkerMemoryBytes(numScenarios, ...
    numBixels,batchSize,cfg)
bytesPerDouble = 8;
kMax = min([cfg.KMax,numBixels,numScenarios]);

scenarioRowsBytes = numScenarios * batchSize * numBixels * bytesPerDouble;
centerBlockBytes = batchSize * numBixels * bytesPerDouble;
perVoxelWorkspaceBytes = estimateOARRadiusFactorMemoryMB(numScenarios,numBixels,cfg) * 1e6;
batchFactorBytes = batchSize * (3*numBixels*kMax + kMax^2) * bytesPerDouble;

workerMemoryBytes = scenarioRowsBytes + centerBlockBytes + ...
    perVoxelWorkspaceBytes + batchFactorBytes;
end

function broadcastBytes = estimateOARRadiusFactorParallelBroadcastBytes(quantity,scenarioMaps)
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

function guardOARRadiusFactorMemory(numScenarios,numBixels,cfg,matRad_cfg)
memoryLimitMB = resolveMemoryLimitMB(cfg);
estimatedMB = estimateOARRadiusFactorMemoryMB(numScenarios,numBixels,cfg);

matRad_cfg.dispInfo(['matRad: Estimated INTERVAL3 OAR radius factor memory ', ...
    'per voxel is %.2f MB (limit %.2f MB).\n'],estimatedMB,memoryLimitMB);

if estimatedMB > memoryLimitMB
    matRad_cfg.dispError(['INTERVAL3 OAR radius factor estimated memory per ', ...
        'voxel is %.2f MB, which exceeds MemoryLimitMB %.2f MB. Increase ', ...
        'cfg.MemoryLimitMB, reduce cfg.KMax, reduce the number of scenarios ', ...
        'or bixels, or use INTERVAL2.'],estimatedMB,memoryLimitMB);
end
end

function estimatedMB = estimateOARRadiusFactorMemoryMB(numScenarios,numBixels,cfg)
bytesPerDouble = 8;
kMax = min([cfg.KMax,numBixels,numScenarios]);

scenarioMatrixBytes = 3 * numScenarios * numBixels * bytesPerDouble;
scenarioGramBytes = 3 * numScenarios^2 * bytesPerDouble;
factorBytes = (3*numBixels*kMax + kMax^2) * bytesPerDouble;

estimatedMB = (scenarioMatrixBytes + scenarioGramBytes + ...
    factorBytes) / 1e6;
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

function [factor,rank] = buildWeightedScenarioRadiusFactor(centeredScenarioMatrix, ...
    scenarioWeights,cfg,numBixels)
if strcmp(cfg.RadiusMode,'extreme')
    [factor,rank] = buildExtremeScenarioRadiusFactor( ...
        centeredScenarioMatrix,scenarioWeights,numBixels);
else
    [factor,rank] = buildStdScenarioRadiusFactor(centeredScenarioMatrix, ...
        scenarioWeights,cfg,numBixels);
end
end

function [factor,rank] = buildStdScenarioRadiusFactor(centeredScenarioMatrix, ...
    scenarioWeights,cfg,numBixels)
tolerance = 1e-12;

kMax = min(cfg.KMax,numBixels);
if kMax == 0
    [factor,rank] = zeroRankRadiusFactor(numBixels);
    return;
end

positiveWeightIx = scenarioWeights(:) > 0;
scenarioWeights = scenarioWeights(positiveWeightIx);
centeredScenarioMatrix = centeredScenarioMatrix(positiveWeightIx,:);

activeBixels = find(any(centeredScenarioMatrix,1));
if isempty(activeBixels)
    [factor,rank] = zeroRankRadiusFactor(numBixels);
    return;
end

centeredScenarioMatrix = centeredScenarioMatrix(:,activeBixels);
numScenarios = size(centeredScenarioMatrix,1);
weightedScenarioMatrix = spdiags(sqrt(scenarioWeights),0, ...
    numScenarios,numScenarios) * centeredScenarioMatrix;

if nnz(weightedScenarioMatrix) == 0 || ...
        norm(weightedScenarioMatrix,'fro') <= tolerance
    [factor,rank] = zeroRankRadiusFactor(numBixels);
    return;
end

% X'*X and X*X' share non-zero eigenvalues; diagonalize the smaller
% scenario Gram and expand compact bixel-space factors.
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
    [factor,rank] = zeroRankRadiusFactor(numBixels);
    return;
end

eigenValues = eigenValues(positiveIx);
leftEigenVectors = leftEigenVectors(:,positiveIx);

if strcmp(cfg.KMode,'dynamic')
    rank = selectScenarioFactorRank(eigenValues,kMax, ...
        cfg.RetentionThreshold,tolerance);
else
    rank = min(kMax,numel(eigenValues));
end

if rank == 0
    [factor,rank] = zeroRankRadiusFactor(numBixels);
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

function [factor,rank] = buildExtremeScenarioRadiusFactor(centeredScenarioMatrix, ...
    scenarioWeights,numBixels)
tolerance = 1e-12;

positiveWeightIx = scenarioWeights(:) > 0;
centeredScenarioMatrix = centeredScenarioMatrix(positiveWeightIx,:);
if isempty(centeredScenarioMatrix)
    [factor,rank] = zeroRankRadiusFactor(numBixels);
    return;
end

delta = max(abs(centeredScenarioMatrix),[],1);
delta(abs(delta) <= tolerance) = 0;
activeBixels = find(delta);

rank = numel(activeBixels);
if rank == 0
    [factor,rank] = zeroRankRadiusFactor(numBixels);
    return;
end

deltaValues = full(delta(activeBixels(:)));
factor = sparse(activeBixels(:),(1:rank)',deltaValues(:), ...
    double(numBixels),double(rank));
end

function rank = selectScenarioFactorRank(eigenValues,kMax,retentionThreshold,tolerance)
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

function [factor,rank] = zeroRankRadiusFactor(numBixels)
rank = 0;
factor = sparse(numBixels,0);
end
