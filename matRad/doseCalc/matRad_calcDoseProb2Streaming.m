function [pln_prob2,dij_prob2Context] = matRad_calcDoseProb2Streaming(ct,cst,stf,pln,varargin)
% matRad_calcDoseProb2Streaming calculates PROB2 data scenario by scenario
%
% call
%   [pln_prob2,dij_prob2Context] = matRad_calcDoseProb2Streaming(ct,cst,stf,pln)
%   [pln_prob2,dij_prob2Context] = matRad_calcDoseProb2Streaming(ct,cst,stf,pln,dij)
%   [pln_prob2,dij_prob2Context] = matRad_calcDoseProb2Streaming(ct,cst,stf,pln,cfg)
%   [pln_prob2,dij_prob2Context] = matRad_calcDoseProb2Streaming(ct,cst,stf,pln,dij,cfg)
%
% input
%   ct:     matRad ct struct
%   cst:    matRad cst cell array
%   stf:    matRad steering information struct
%   pln:    matRad pln struct with robust scenario model
%   dij:    optional precomputed robust dose influence struct. If provided,
%           streaming uses dij scenario matrices instead of recalculating
%           scenario dose influence data
%   cfg:    configuration struct. Optional fields accepted by
%           matRad_calcDoseProb2 plus:
%           UseParallel: use safe available scenario parallelism for the
%               first pass and Omega pass when the Parallel Computing
%               Toolbox and enough workers/memory are available.
%               matRad may create or reduce the active parallel pool.
%               Precomputed dij inputs fall back to serial streaming.
%           SecondPassStrategy: 'disk' (default) or 'recompute'
%           CacheRoot: root folder for disk blocks
%           KeepCache: keep the hash cache folder after the run (default false)
%
% output
%   pln_prob2:        plan struct containing selected-VOI propOpt.dij_prob2.
%                     The propOpt.dij_prob2.streamingSize field summarizes
%                     compact result bytes and peak streaming auxiliary
%                     bytes used during precomputation.
%   dij_prob2Context: lightweight single-scenario dij context for
%                     probabilistic fluence optimization
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
[~,cfg] = matRad_parseScenarioDoseStreamingArguments(matRad_cfg,varargin{:});
timer = tic;
cacheContext = [];

try
    cfg = normalizeProb2StreamingConfig(cfg,matRad_cfg);
    scenarioInfo = matRad_resolveStreamingScenarioSelection(ct,pln,cfg,matRad_cfg);
    [provider,dijForResolve] = matRad_initializeScenarioRowProvider( ...
        ct,cst,stf,pln,cfg,scenarioInfo,matRad_cfg,'PROB2');
    provider.sizeTelemetry = matRad_initializeScenarioDoseSizeTelemetry();

    ctx = matRad_resolveScenarioDoseInputs(ct,cst,pln,dijForResolve,cfg, ...
        'PROB2',matRad_cfg);
    ctx.scenarioIds = scenarioInfo.scenarioIds(:);
    if isfield(dijForResolve,'doseGrid')
        ctx.doseGrid = dijForResolve.doseGrid;
    end
    if isfield(dijForResolve,'ctGrid')
        ctx.ctGrid = dijForResolve.ctGrid;
    end

    cfg = ctx.cfg;
    cfg.SecondPassStrategy = provider.secondPassStrategy;
    cfg.KeepCache = provider.keepCache;
    cfg.CacheRoot = provider.cacheRoot;
    quantity = validateProb2StreamingQuantity(ctx.quantity,matRad_cfg);

    matRad_cfg.dispInfo(['matRad: Calculating PROB2 streaming data for ', ...
        'quantity ''%s'' and %d scenarios.\n'], ...
        quantity.name,numel(ctx.scenarioDijIx));
    matRad_cfg.dispInfo(['matRad: PROB2 streaming reference CT scenario ', ...
        '%d, %d voxels, %d bixels.\n'], ...
        cfg.refScen,ctx.numVoxels,ctx.numBixels);

    expectedRows = resolveProb2StreamingExpectedRows(ctx);
    expectedBatches = matRad_makeScenarioDoseRowBatches(expectedRows, ...
        matRad_resolveScenarioDoseBatchSize(numel(expectedRows), ...
        numel(ctx.scenarioDijIx),ctx.numBixels,cfg));
    voiRows = resolveProb2StreamingVoiRows(ctx);
    voiBatches = makeProb2StreamingVoiBatches(voiRows,ctx,cfg);

    needsSecondPass = hasProb2StreamingVoiRows(voiRows);
    needsDiskCache = strcmp(cfg.SecondPassStrategy,'disk') && needsSecondPass;
    if needsDiskCache
        signatureExtras = struct();
        signatureExtras.voiRows = voiRows;
        cacheContext = matRad_initializeScenarioDoseCache(cfg,ctx,quantity, ...
            stf,matRad_cfg,'PROB2',signatureExtras);
        provider.cacheContext = cacheContext;
    end

    firstPassWork = struct();
    firstPassWork.primaryBatches = expectedBatches;
    firstPassWork.secondaryBatches = voiBatches;
    firstPassWork.cachePrimaryRows = false;
    firstPassWork.cacheSecondaryRows = needsDiskCache;

    dij_prob2 = initializeProb2StreamingStruct(ctx,quantity);
    dij_prob2.voiSubIx = voiRows;
    dij_prob2.precomputeMode = 'streaming';
    dij_prob2.secondPassStrategy = cfg.SecondPassStrategy;

    [dij_prob2,provider] = accumulateProb2StreamingFirstPass( ...
        dij_prob2,provider,ctx,quantity,cfg,firstPassWork,matRad_cfg);

    if needsSecondPass
        [dij_prob2.Omega,provider,omegaParallel] = accumulateProb2StreamingOmega( ...
            dij_prob2,provider,ctx,quantity,cfg,voiRows,voiBatches, ...
            matRad_cfg);
        if cfg.CollectTiming
            dij_prob2.timing.parallelScenario.omega = omegaParallel;
        end
    else
        dij_prob2.Omega = cell(size(voiRows));
    end

    if cfg.CollectTiming
        dij_prob2.timing.totalSeconds = toc(timer);
    end

    if needsDiskCache && cfg.KeepCache
        dij_prob2.cacheDir = cacheContext.runDir;
    end
    dij_prob2.streamingSize = matRad_buildScenarioDoseStreamingSize( ...
        dij_prob2,provider,cfg);

    pln_prob2 = pln;
    if ~isfield(pln_prob2,'propOpt') || ~isstruct(pln_prob2.propOpt)
        pln_prob2.propOpt = struct();
    end
    pln_prob2.propOpt.dij_prob2 = dij_prob2;
    dij_prob2Context = buildProb2StreamingDijContext(ct,cst,stf,pln, ...
        dij_prob2,quantity,cfg);
    pln_prob2.multScen = dij_prob2Context.scenarioModel;

    matRad_cleanupScenarioDoseCache(cacheContext,cfg.KeepCache);
    matRad_cfg.dispInfo('matRad: Finished PROB2 streaming calculation in %.2f s.\n', ...
        toc(timer));
catch ME
    if ~isempty(cacheContext)
        matRad_cleanupScenarioDoseCache(cacheContext,false);
    end
    rethrow(ME);
end

end

function cfg = normalizeProb2StreamingConfig(cfg,matRad_cfg)
cfg = matRad_normalizeScenarioDoseStreamingConfig(cfg,matRad_cfg, ...
    fullfile(matRad_cfg.matRadRoot,'cache','doseProb2'));
end

function quantity = validateProb2StreamingQuantity(quantity,matRad_cfg)
if strcmpi(quantity.name,'physicalDose')
    return;
end

if strcmpi(quantity.name,'RBExD') && strcmp(quantity.field,'physicalDose') && ...
        quantity.scale ~= 1
    return;
end

matRad_cfg.dispError(['PROB2 V1 supports only physicalDose and RBExD ', ...
    'with scalar constant dij.RBE.']);
end

function dij_prob2 = initializeProb2StreamingStruct(ctx,quantity)
dij_prob2 = struct();
dij_prob2.expected = sparse(ctx.numVoxels,ctx.numBixels);
dij_prob2.Omega = cell(0,1);
dij_prob2.voiSubIx = cell(0,1);
dij_prob2.quantity = quantity.name;
dij_prob2.quantityField = quantity.field;
dij_prob2.quantityScale = quantity.scale;
dij_prob2.optimizationQuantity = quantity.optimizationQuantity;
dij_prob2.refScen = ctx.cfg.refScen;
dij_prob2.scenarioDijIx = ctx.scenarioDijIx(:);
dij_prob2.scenarioCtScenIds = ctx.scenarioCtScenIds(:);
dij_prob2.scenarioWeights = ctx.scenarioWeights(:);
dij_prob2.probabilisticMode = 'PROB2';
if isfield(ctx.cfg,'CollectTiming') && ctx.cfg.CollectTiming
    dij_prob2.timing = struct();
    dij_prob2.timing.useParallelRequested = logical(ctx.cfg.UseParallel);
    dij_prob2.timing.numScenarios = numel(ctx.scenarioDijIx);
    dij_prob2.timing.numBixels = ctx.numBixels;
    dij_prob2.timing.totalSeconds = 0;
    dij_prob2.timing.parallelScenario = struct('firstPass',false, ...
        'omega',false);
end
end

function rows = resolveProb2StreamingExpectedRows(ctx)
rows = unique([ctx.targetRows(:); ctx.oarRows(:)],'stable');
end

function voiRows = resolveProb2StreamingVoiRows(ctx)
voiRows = cell(size(ctx.cstDoseGrid,1),1);
selectedStructures = [ctx.targetStructures(:); ctx.oarStructures(:)];

for i = 1:numel(selectedStructures)
    rowIx = selectedStructures(i).cstRow;
    voxels = selectedStructures(i).voxelIx(:);
    if isempty(voiRows{rowIx})
        voiRows{rowIx} = voxels;
    else
        voiRows{rowIx} = unique([voiRows{rowIx}; voxels],'stable');
    end
end
end

function dij_prob2Context = buildProb2StreamingDijContext(ct,cst,stf,pln, ...
    dij_prob2,quantity,cfg)
numBixels = size(dij_prob2.expected,2);
robustDij = [];
if isfield(cfg,'PrecomputedDij')
    robustDij = cfg.PrecomputedDij;
end
dij_prob2Context = matRad_buildNominalDijContext( ...
    ct,cst,stf,pln,cfg,numBixels,robustDij);
if strcmpi(quantity.name,'RBExD') && isfield(dij_prob2Context,'RBE')
    dij_prob2Context.RBE = 1;
end
dij_prob2Context.probabilisticQuantity = quantity.name;
dij_prob2Context.probabilisticQuantityField = quantity.field;
end

function [dij_prob2,provider] = accumulateProb2StreamingFirstPass( ...
    dij_prob2,provider,ctx,quantity,cfg,firstPassWork,matRad_cfg)
numScenarios = numel(ctx.scenarioDijIx);
expectedBatches = firstPassWork.primaryBatches;
voiBatches = firstPassWork.secondaryBatches;
cacheVoiRows = firstPassWork.cacheSecondaryRows;
numExpectedRows = sum(cellfun(@(batch) numel(batch.rows),expectedBatches));

matRad_cfg.dispInfo(['matRad: Streaming PROB2 first pass over %d ', ...
    'scenario(s) and %d selected dose voxels.\n'], ...
    numScenarios,numExpectedRows);

[useParallel,parallelProvider] = matRad_configureScenarioDoseStreamingParallel( ...
    provider,ctx,cfg,matRad_cfg,'streaming PROB2 first pass',[]);
if useParallel
    expectedRows = sparse(numExpectedRows,ctx.numBixels);
    logLevel = matRad_cfg.logLevel;
    progressReporter = matRad_createScenarioDoseProgressReporter(matRad_cfg, ...
        cfg,'Expected influence',numScenarios*numel(expectedBatches),true);
    progressCleanup = onCleanup(progressReporter.cleanup);
    progressQueue = progressReporter.queue;
    scenarioResults = cell(numScenarios,1);
    parfor s = 1:numScenarios
        workerCfg = MatRad_Config.instance();
        workerCfg.logLevel = logLevel;
        scenarioResults{s} = computeProb2StreamingFirstPassScenario( ...
            parallelProvider,ctx,quantity,cfg,expectedBatches,voiBatches, ...
            cacheVoiRows,s,numExpectedRows,workerCfg,progressQueue, ...
            numScenarios);
    end

    telemetryBlocks = cell(numScenarios,1);
    for s = 1:numScenarios
        result = scenarioResults{s};
        expectedRows = expectedRows + result.expectedRows;
        telemetryBlocks{s} = result.sizeTelemetry;
    end
    provider = matRad_mergeScenarioDoseSizeTelemetry(provider,telemetryBlocks);
    if cfg.CollectTiming
        dij_prob2.timing.parallelScenario.firstPass = true;
    end

    for b = 1:numel(expectedBatches)
        batch = expectedBatches{b};
        dij_prob2.expected(batch.rows,:) = expectedRows(batch.localIx,:);
    end
    return;
end

progressReporter = matRad_createScenarioDoseProgressReporter(matRad_cfg, ...
    cfg,'Expected influence',numScenarios*numel(expectedBatches),false);
for s = 1:numScenarios
    [provider,source] = matRad_beginScenarioDoseRowsProvider(provider,ctx, ...
        quantity,s,matRad_cfg);
    scenarioWeight = ctx.scenarioWeights(s);

    for b = 1:numel(expectedBatches)
        batch = expectedBatches{b};
        rows = matRad_getScenarioDoseProviderRows(source,ctx.scenarioMaps{s}, ...
            batch.rows,matRad_cfg);
        dij_prob2.expected(batch.rows,:) = dij_prob2.expected(batch.rows,:) + ...
            scenarioWeight .* rows;
        progressReporter.update(sprintf('scenario %d/%d, batch %d/%d', ...
            s,numScenarios,b,numel(expectedBatches)));
    end

    if cacheVoiRows
        for voiIx = 1:numel(voiBatches)
            batches = voiBatches{voiIx};
            for b = 1:numel(batches)
                batch = batches{b};
                rows = matRad_getScenarioDoseProviderRows(source, ...
                    ctx.scenarioMaps{s},batch.rows,matRad_cfg);
                blockBytes = matRad_writeScenarioDoseCacheBlock(provider.cacheContext,s, ...
                    prob2StreamingVoiCacheKind(voiIx),b,batch.rows,rows);
                provider = matRad_updateScenarioDoseDiskCachePeak(provider, ...
                    blockBytes);
            end
        end
    end

    [provider,source] = matRad_endScenarioDoseRowsProvider(provider,source);
end
end

function result = computeProb2StreamingFirstPassScenario( ...
    provider,ctx,quantity,cfg,expectedBatches,voiBatches,cacheVoiRows, ...
    scenarioIx,numExpectedRows,matRad_cfg,progressQueue,numScenarios)
localProvider = provider;
if ~isfield(localProvider,'sizeTelemetry') || isempty(localProvider.sizeTelemetry)
    localProvider.sizeTelemetry = matRad_initializeScenarioDoseSizeTelemetry();
end

expectedRows = sparse(numExpectedRows,ctx.numBixels);
[localProvider,source] = matRad_beginScenarioDoseRowsProvider( ...
    localProvider,ctx,quantity,scenarioIx,matRad_cfg);
scenarioWeight = ctx.scenarioWeights(scenarioIx);

for b = 1:numel(expectedBatches)
    batch = expectedBatches{b};
    rows = matRad_getScenarioDoseProviderRows(source, ...
        ctx.scenarioMaps{scenarioIx},batch.rows,matRad_cfg);
    expectedRows(batch.localIx,:) = expectedRows(batch.localIx,:) + ...
        scenarioWeight .* rows;
    if ~isempty(progressQueue)
        send(progressQueue,sprintf('scenario %d/%d, batch %d/%d', ...
            scenarioIx,numScenarios,b,numel(expectedBatches)));
    end
end

if cacheVoiRows
    for voiIx = 1:numel(voiBatches)
        batches = voiBatches{voiIx};
        for b = 1:numel(batches)
            batch = batches{b};
            rows = matRad_getScenarioDoseProviderRows(source, ...
                ctx.scenarioMaps{scenarioIx},batch.rows,matRad_cfg);
            blockBytes = matRad_writeScenarioDoseCacheBlock( ...
                localProvider.cacheContext,scenarioIx, ...
                prob2StreamingVoiCacheKind(voiIx),b,batch.rows,rows);
            localProvider = matRad_updateScenarioDoseDiskCachePeak( ...
                localProvider,blockBytes);
        end
    end
end

[localProvider,source] = matRad_endScenarioDoseRowsProvider(localProvider,source);

result = struct();
result.expectedRows = expectedRows;
result.sizeTelemetry = localProvider.sizeTelemetry;
end

function [Omega,provider,parallelEnabled] = accumulateProb2StreamingOmega( ...
    dij_prob2,provider,ctx,quantity,cfg,voiRows,voiBatches,matRad_cfg)
Omega = cell(size(voiRows));
numScenarios = numel(ctx.scenarioDijIx);
numBixels = ctx.numBixels;
parallelEnabled = false;

[useParallel,parallelProvider] = matRad_configureScenarioDoseStreamingParallel( ...
    provider,ctx,cfg,matRad_cfg,'streaming PROB2 Omega',[]);
if useParallel
    logLevel = matRad_cfg.logLevel;
    numOmegaBatches = sum(cellfun(@numel,voiBatches));
    progressReporter = matRad_createScenarioDoseProgressReporter(matRad_cfg, ...
        cfg,'Omega',numScenarios*numOmegaBatches,true);
    progressCleanup = onCleanup(progressReporter.cleanup);
    progressQueue = progressReporter.queue;
    scenarioResults = cell(numScenarios,1);
    parfor s = 1:numScenarios
        workerCfg = MatRad_Config.instance();
        workerCfg.logLevel = logLevel;
        scenarioResults{s} = computeProb2StreamingOmegaScenario( ...
            parallelProvider,dij_prob2,ctx,quantity,cfg,voiBatches,s, ...
            workerCfg,progressQueue,numScenarios);
    end

    telemetryBlocks = cell(numScenarios,1);
    for voiIx = 1:numel(voiRows)
        if isempty(voiRows{voiIx})
            Omega{voiIx} = [];
        else
            Omega{voiIx} = sparse(numBixels,numBixels);
        end
    end
    for s = 1:numScenarios
        result = scenarioResults{s};
        for voiIx = 1:numel(voiRows)
            if ~isempty(result.Omega{voiIx})
                Omega{voiIx} = Omega{voiIx} + result.Omega{voiIx};
            end
        end
        telemetryBlocks{s} = result.sizeTelemetry;
    end
    provider = matRad_mergeScenarioDoseSizeTelemetry(provider,telemetryBlocks);
    for voiIx = 1:numel(Omega)
        if ~isempty(Omega{voiIx})
            Omega{voiIx} = sparse(0.5.*(Omega{voiIx} + Omega{voiIx}'));
        end
    end
    parallelEnabled = true;
    return;
end

numOmegaBatches = sum(cellfun(@numel,voiBatches));
progressReporter = matRad_createScenarioDoseProgressReporter(matRad_cfg, ...
    cfg,'Omega',numScenarios*numOmegaBatches,false);
for voiIx = 1:numel(voiRows)
    if isempty(voiRows{voiIx})
        Omega{voiIx} = [];
        continue;
    end

    matRad_cfg.dispInfo(['matRad: Streaming PROB2 second pass Omega for ', ...
        'VOI %d (%d voxels).\n'],voiIx,numel(voiRows{voiIx}));
    omega = sparse(numBixels,numBixels);
    batches = voiBatches{voiIx};

    if strcmp(cfg.SecondPassStrategy,'disk')
        for b = 1:numel(batches)
            batch = batches{b};
            centeredRows = cell(numScenarios,1);
            for s = 1:numScenarios
                rows = matRad_readScenarioDoseCacheBlock(provider.cacheContext, ...
                    s,prob2StreamingVoiCacheKind(voiIx),b,batch.rows);
                centeredRows{s} = rows - dij_prob2.expected(batch.rows,:);
                progressReporter.update(sprintf( ...
                    'scenario %d/%d, VOI %d, batch %d/%d', ...
                    s,numScenarios,voiIx,b,numel(batches)));
            end
            omega = omega + accumulateProb2StreamingCenteredOmegaBatch(centeredRows, ...
                ctx.scenarioWeights,numBixels);
        end
    else
        for s = 1:numScenarios
            [provider,source] = matRad_beginScenarioDoseRowsProvider(provider, ...
                ctx,quantity,s,matRad_cfg);
            for b = 1:numel(batches)
                batch = batches{b};
                rows = matRad_getScenarioDoseProviderRows(source, ...
                    ctx.scenarioMaps{s},batch.rows,matRad_cfg);
                centeredRows = {rows - dij_prob2.expected(batch.rows,:)};
                provider = matRad_updateScenarioDoseMemoryPeak(provider, ...
                    source,rows,centeredRows);
                omega = omega + accumulateProb2StreamingCenteredOmegaBatch(centeredRows, ...
                    ctx.scenarioWeights(s),numBixels);
                progressReporter.update(sprintf( ...
                    'scenario %d/%d, VOI %d, batch %d/%d', ...
                    s,numScenarios,voiIx,b,numel(batches)));
            end
            [provider,source] = matRad_endScenarioDoseRowsProvider(provider,source);
        end
    end

    Omega{voiIx} = sparse(0.5.*(omega + omega'));
end
end

function result = computeProb2StreamingOmegaScenario( ...
    provider,dij_prob2,ctx,quantity,cfg,voiBatches,scenarioIx,matRad_cfg, ...
    progressQueue,numScenarios)
localProvider = provider;
if ~isfield(localProvider,'sizeTelemetry') || isempty(localProvider.sizeTelemetry)
    localProvider.sizeTelemetry = matRad_initializeScenarioDoseSizeTelemetry();
end

Omega = cell(size(dij_prob2.voiSubIx));
numBixels = ctx.numBixels;
source = [];
if strcmp(cfg.SecondPassStrategy,'recompute')
    [localProvider,source] = matRad_beginScenarioDoseRowsProvider( ...
        localProvider,ctx,quantity,scenarioIx,matRad_cfg);
end

for voiIx = 1:numel(voiBatches)
    batches = voiBatches{voiIx};
    if isempty(batches)
        Omega{voiIx} = [];
        continue;
    end

    omega = sparse(numBixels,numBixels);
    for b = 1:numel(batches)
        batch = batches{b};
        if strcmp(cfg.SecondPassStrategy,'disk')
            rows = matRad_readScenarioDoseCacheBlock(localProvider.cacheContext, ...
                scenarioIx,prob2StreamingVoiCacheKind(voiIx),b,batch.rows);
        else
            rows = matRad_getScenarioDoseProviderRows(source, ...
                ctx.scenarioMaps{scenarioIx},batch.rows,matRad_cfg);
        end
        centeredRows = {rows - dij_prob2.expected(batch.rows,:)};
        if strcmp(cfg.SecondPassStrategy,'recompute')
            localProvider = matRad_updateScenarioDoseMemoryPeak( ...
                localProvider,source,rows,centeredRows);
        end
        omega = omega + accumulateProb2StreamingCenteredOmegaBatch( ...
            centeredRows,ctx.scenarioWeights(scenarioIx),numBixels);
        if ~isempty(progressQueue)
            send(progressQueue,sprintf( ...
                'scenario %d/%d, VOI %d, batch %d/%d', ...
                scenarioIx,numScenarios,voiIx,b,numel(batches)));
        end
    end
    Omega{voiIx} = omega;
end

if strcmp(cfg.SecondPassStrategy,'recompute')
    [localProvider,source] = matRad_endScenarioDoseRowsProvider(localProvider,source);
end

result = struct();
result.Omega = Omega;
result.sizeTelemetry = localProvider.sizeTelemetry;
end

function omega = accumulateProb2StreamingCenteredOmegaBatch(centeredRows,scenarioWeights, ...
    numBixels)
omega = sparse(numBixels,numBixels);

for s = 1:numel(centeredRows)
    centeredScenarioRows = centeredRows{s};
    activeBixels = find(any(centeredScenarioRows,1));
    if isempty(activeBixels)
        continue;
    end

    centeredActive = centeredScenarioRows(:,activeBixels);
    omega(activeBixels,activeBixels) = ...
        omega(activeBixels,activeBixels) + ...
        centeredActive' * (scenarioWeights(s).*centeredActive);
end
end

function voiBatches = makeProb2StreamingVoiBatches(voiRows,ctx,cfg)
voiBatches = cell(size(voiRows));
for voiIx = 1:numel(voiRows)
    rows = voiRows{voiIx};
    if isempty(rows)
        voiBatches{voiIx} = {};
        continue;
    end
    batchSize = matRad_resolveScenarioDoseBatchSize(numel(rows), ...
        numel(ctx.scenarioDijIx),ctx.numBixels,cfg);
    voiBatches{voiIx} = matRad_makeScenarioDoseRowBatches(rows,batchSize);
end
end

function tf = hasProb2StreamingVoiRows(values)
tf = any(~cellfun(@isempty,values(:)));
end

function kind = prob2StreamingVoiCacheKind(voiIx)
kind = sprintf('voi_%04d',voiIx);
end
