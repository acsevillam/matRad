function [plnInterval, dijIntervalContext] = matRad_calcDoseInterval(ct, cst, stf, pln, varargin)
% matRad_calcDoseInterval calculates interval data scenario by scenario
%
% call
%   [plnInterval,dijIntervalContext] = matRad_calcDoseInterval(ct,cst,stf,pln,cfg)
%   [plnInterval,dijIntervalContext] = matRad_calcDoseInterval(ct,cst,stf,pln,dij,cfg)
%
% input
%   ct:     matRad ct struct
%   cst:    matRad cst cell array
%   stf:    matRad steering information struct
%   pln:    matRad pln struct with robust scenario model
%   dij:    optional precomputed robust dose influence struct. If provided,
%           precomputed scenario matrices from dij are used instead of
%           recalculating scenario dose influence data.
%   cfg:    configuration struct. Required field:
%           IntervalMode: 'INTERVAL2' or 'INTERVAL3'
%           Optional fields:
%           UseParallel: use safe available scenario parallelism for
%               scenario-reducible passes when the Parallel
%               Computing Toolbox and enough workers/memory are available.
%               If target extreme radius parallel aggregation exceeds
%               MemoryLimitMB but serial target extreme radius fits, only
%               that second pass falls back to serial execution.
%               matRad may create, reduce, or increase the active parallel
%               pool.
%               INTERVAL3 OAR radius factors keep their batch/voxel path.
%               Precomputed dij inputs run serially.
%           SecondPassStrategy: 'disk' (default) or 'recompute'
%           CacheRoot: root folder for disk blocks. Defaults to a temporary
%               folder outside the matRad checkout.
%           KeepCache: keep the hash cache folder after the run (default false)
%           MemoryLimitMB: positive scalar MB or 'auto' (default). Auto uses
%               scheduler/cgroup memory first and falls back conservatively.
%           MemoryLimitFraction: fraction of detected memory used by 'auto'
%           MemoryLimitFallbackMB: fallback MB when 'auto' cannot detect memory
%
% output
%   plnInterval:        plan struct containing propOpt.dij_interval. The
%                        propOpt.dij_interval.precomputeSize field summarizes
%                        compact result bytes and peak auxiliary bytes used
%                        during precomputation.
%   dijIntervalContext: lightweight single-scenario dij context for
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

matRadCfg = MatRad_Config.instance();
[~, cfg] = ScenarioBatch.Config.matRad_parseScenarioDosePrecomputeArguments(matRadCfg, varargin{:});
if ~isfield(cfg, 'IntervalMode') || isempty(cfg.IntervalMode)
    matRadCfg.dispError('cfg.IntervalMode must be ''INTERVAL2'' or ''INTERVAL3''.');
end
intervalMode = cfg.IntervalMode;

timer = tic;
cacheContext = [];

try
    [cfg, intervalMode] = DoseInterval.matRad_normalizeDoseIntervalConfig(cfg, intervalMode, matRadCfg);
    scenarioInfo = ScenarioBatch.Input.matRad_selectScenarioDoseRealizations(ct, pln, cfg, matRadCfg);
    [provider, dijForResolve] = ScenarioBatch.Source.matRad_createScenarioDoseRowProvider( ...
                                                                                          ct, cst, stf, pln, cfg, ...
                                                                                          scenarioInfo, matRadCfg, 'dose interval');
    provider.resourceUsage = ScenarioBatch.Resources.matRad_buildScenarioDoseResourceUsage();

    ctx = ScenarioBatch.Input.matRad_buildScenarioDoseContext(ct, cst, pln, dijForResolve, cfg, ...
                                                              intervalMode, matRadCfg, scenarioInfo);
    ctx.scenarioIds = scenarioInfo.scenarioIds(:);
    if isfield(dijForResolve, 'doseGrid')
        ctx.doseGrid = dijForResolve.doseGrid;
    end
    if isfield(dijForResolve, 'ctGrid')
        ctx.ctGrid = dijForResolve.ctGrid;
    end
    cfg = ctx.cfg;
    cfg.IntervalMode = intervalMode;
    cfg.SecondPassStrategy = provider.secondPassStrategy;
    cfg.KeepCache = provider.keepCache;
    cfg.CacheRoot = provider.cacheRoot;
    quantity = ctx.quantity;

    matRadCfg.dispInfo(['matRad: Calculating %s dose interval using ', ...
                        'scenario-batch precomputation for quantity ''%s'' and %d scenarios.\n'], ...
                       intervalMode, quantity.name, numel(ctx.scenarioDijIx));

    [targetBatches, oarBatches] = DoseInterval.matRad_buildDoseIntervalBatches(ctx, cfg);

    needsTargetSecondPass = ~isempty(ctx.targetRows) && ...
        strcmp(cfg.RadiusMode, 'extreme');
    needsOARSecondPass = ~isempty(ctx.oarRows) && strcmp(intervalMode, 'INTERVAL3');
    needsDiskCache = strcmp(cfg.SecondPassStrategy, 'disk') && ...
        (needsTargetSecondPass || needsOARSecondPass);

    if needsOARSecondPass
        DoseInterval.matRad_guardDoseIntervalOARRadiusFactorMemory(numel(ctx.scenarioDijIx), ...
                                                                   ctx.numBixels, cfg, matRadCfg);
    end

    if needsDiskCache
        signatureExtras = struct();
        signatureExtras.intervalMode = cfg.IntervalMode;
        signatureExtras.radiusMode = cfg.RadiusMode;
        signatureExtras.targetRows = ctx.targetRows(:);
        signatureExtras.oarRows = ctx.oarRows(:);
        cacheContext = ScenarioBatch.Cache.matRad_createScenarioDoseCache(cfg, ctx, quantity, ...
                                                                          stf, matRadCfg, 'dose interval', signatureExtras);
        provider.cacheContext = cacheContext;
    end

    firstPass = struct();
    firstPass.targetBatches = targetBatches;
    firstPass.oarBatches = oarBatches;
    firstPass.cacheTargetRows = needsTargetSecondPass && ...
        strcmp(cfg.SecondPassStrategy, 'disk');
    firstPass.cacheOARRows = needsOARSecondPass && ...
        strcmp(cfg.SecondPassStrategy, 'disk');

    dijInterval = DoseInterval.matRad_buildDoseIntervalDij(ctx, quantity, cfg, intervalMode);
    [dijInterval, provider] = matRad_accumulateIntervalScenarioBatchFirstPass( ...
                                                                              dijInterval, provider, ctx, quantity, cfg, firstPass, matRadCfg);

    if needsTargetSecondPass
        [dijInterval, provider] = matRad_accumulateIntervalScenarioBatchTargetExtremeRadius( ...
                                                                                            dijInterval, provider, ctx, quantity, ...
                                                                                            cfg, targetBatches, matRadCfg);
    end

    if needsOARSecondPass
        [dijInterval, provider] = matRad_accumulateIntervalScenarioBatchOARRadiusFactors( ...
                                                                                         dijInterval, provider, ctx, quantity, ...
                                                                                         cfg, oarBatches, matRadCfg);
    end

    if cfg.CollectTiming
        dijInterval.timing.totalSeconds = toc(timer);
    end

    if needsDiskCache && cfg.KeepCache
        dijInterval.cacheDir = cacheContext.runDir;
    end
    dijInterval.precomputeSize = ScenarioBatch.Resources.matRad_buildScenarioDosePrecomputeSize( ...
                                                                                                dijInterval, provider, cfg);

    plnInterval = pln;
    if ~isfield(plnInterval, 'propOpt') || ~isstruct(plnInterval.propOpt)
        plnInterval.propOpt = struct();
    end
    plnInterval.propOpt.dij_interval = dijInterval;
    dijIntervalContext = DoseInterval.matRad_buildDoseIntervalDijContext(ct, cst, stf, pln, ...
                                                                         dijInterval, quantity, cfg);
    plnInterval.multScen = dijIntervalContext.scenarioModel;

    ScenarioBatch.Cache.matRad_cleanupScenarioDoseCache(cacheContext, cfg.KeepCache);
    matRadCfg.dispInfo('matRad: Finished %s scenario-batch dose interval calculation in %.2f s.\n', ...
                       intervalMode, toc(timer));
catch ME
    if ~isempty(cacheContext)
        ScenarioBatch.Cache.matRad_cleanupScenarioDoseCache(cacheContext, false);
    end
    rethrow(ME);
end

end

function [dijInterval, provider] = matRad_accumulateIntervalScenarioBatchFirstPass( ...
                                                                                   dijInterval, provider, ctx, quantity, cfg, ...
                                                                                   firstPass, matRadCfg)
numScenarios = numel(ctx.scenarioDijIx);
targetBatches = firstPass.targetBatches;
oarBatches = firstPass.oarBatches;
cacheTargetRows = firstPass.cacheTargetRows;
cacheOARRows = firstPass.cacheOARRows;
targetCenter = sparse(numel(ctx.targetRows), ctx.numBixels);
oarCenter = sparse(numel(ctx.oarRows), ctx.numBixels);
targetSecondMoment = sparse(ctx.numBixels, ctx.numBixels);

matRadCfg.dispInfo(['matRad: Scenario-batch first pass over %d scenario(s), ', ...
                    '%d target voxels, and %d OAR voxels.\n'], numScenarios, ...
                   numel(ctx.targetRows), numel(ctx.oarRows));

stageName = 'scenario-batch interval first pass';
resultBytesPerScenario = matRad_estimateIntervalFirstPassResultBytes(ctx, cfg);
accumulatorBytes = resultBytesPerScenario;
[useParallel, parallelProvider, parallelPlan] = ScenarioBatch.Parallel.matRad_configureScenarioDoseParallel( ...
                                                                                                           provider, ctx, cfg, matRadCfg, ...
                                                                                                           stageName, [], ...
                                                                                                           resultBytesPerScenario, ...
                                                                                                           accumulatorBytes);
if useParallel
    logLevel = matRadCfg.logLevel;
    targetProgressWorkItems = numScenarios * numel(targetBatches);
    targetProgressReporter = ScenarioBatch.Telemetry.matRad_createScenarioDoseProgressReporter( ...
                                                                                               matRadCfg, cfg, 'Target center', ...
                                                                                               targetProgressWorkItems, true);
    targetProgressCleanup = onCleanup(targetProgressReporter.cleanup);
    targetProgressQueue = targetProgressReporter.queue;
    oarProgressWorkItems = numScenarios * numel(oarBatches);
    oarProgressReporter = ScenarioBatch.Telemetry.matRad_createScenarioDoseProgressReporter( ...
                                                                                            matRadCfg, cfg, 'OAR center', ...
                                                                                            oarProgressWorkItems, true);
    oarProgressCleanup = onCleanup(oarProgressReporter.cleanup);
    oarProgressQueue = oarProgressReporter.queue;
    scenarioChunks = ScenarioBatch.Parallel.matRad_buildScenarioChunks( ...
                                                                       numScenarios, parallelPlan.chunkSize);
    for chunkIx = 1:numel(scenarioChunks)
        chunkScenarios = scenarioChunks{chunkIx};
        chunkResults = cell(numel(chunkScenarios), 1);
        parfor chunkLocalIx = 1:numel(chunkScenarios)
            scenarioIx = chunkScenarios(chunkLocalIx);
            workerCfg = MatRad_Config.instance();
            workerCfg.logLevel = logLevel;
            chunkResults{chunkLocalIx} = matRad_computeIntervalScenarioBatchFirstPassScenario( ...
                                                                                              parallelProvider, ctx, quantity, cfg, targetBatches, ...
                                                                                              oarBatches, cacheTargetRows, cacheOARRows, scenarioIx, ...
                                                                                              workerCfg, targetProgressQueue, oarProgressQueue, ...
                                                                                              numScenarios);
        end

        for chunkLocalIx = 1:numel(chunkResults)
            result = chunkResults{chunkLocalIx};
            targetCenter = targetCenter + result.targetCenter;
            oarCenter = oarCenter + result.oarCenter;
            targetSecondMoment = targetSecondMoment + result.targetSecondMoment;
        end
        provider = ScenarioBatch.Resources.matRad_mergeScenarioDoseResourceUsage(provider, chunkResults);
        chunkResults = [];
    end
    if cfg.CollectTiming
        dijInterval.timing.parallelScenario.firstPass = true;
    end

    if ~isempty(ctx.targetRows)
        dijInterval.center(ctx.targetRows, :) = targetCenter;
    end
    if ~isempty(ctx.oarRows)
        dijInterval.center(ctx.oarRows, :) = oarCenter;
    end
    if strcmp(cfg.RadiusMode, 'std')
        dijInterval.radius = targetSecondMoment;
    end
    return
end

targetProgressWorkItems = numScenarios * numel(targetBatches);
targetProgressReporter = ScenarioBatch.Telemetry.matRad_createScenarioDoseProgressReporter( ...
                                                                                           matRadCfg, cfg, 'Target center', ...
                                                                                           targetProgressWorkItems, false);
oarProgressWorkItems = numScenarios * numel(oarBatches);
oarProgressReporter = ScenarioBatch.Telemetry.matRad_createScenarioDoseProgressReporter( ...
                                                                                        matRadCfg, cfg, 'OAR center', ...
                                                                                        oarProgressWorkItems, false);
targetCacheKind = 'target';
oarCacheKind = 'oar';
for s = 1:numScenarios
    [provider, source] = ScenarioBatch.Source.matRad_createScenarioDoseRowSource(provider, ctx, quantity, s, matRadCfg);
    scenarioWeight = ctx.scenarioWeights(s);

    for b = 1:numel(targetBatches)
        batch = targetBatches{b};
        rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ctx.scenarioMaps{s}, ...
                                                               batch.rows, matRadCfg);
        targetCenter(batch.localIx, :) = targetCenter(batch.localIx, :) + ...
            scenarioWeight .* rows;
        if strcmp(cfg.RadiusMode, 'std')
            targetSecondMoment = targetSecondMoment + rows' * ...
                (scenarioWeight .* rows);
        end
        if cacheTargetRows
            blockBytes = ScenarioBatch.Cache.matRad_writeScenarioDoseCacheBlock(provider.cacheContext, s, ...
                                                                                targetCacheKind, b, batch.rows, rows);
            provider = ScenarioBatch.Resources.matRad_updateScenarioDoseDiskCacheUsage(provider, ...
                                                                                       blockBytes);
        end
        targetProgressReporter.update(sprintf( ...
                                              'scenario %d/%d, batch %d/%d', s, numScenarios, b, ...
                                              numel(targetBatches)));
    end

    for b = 1:numel(oarBatches)
        batch = oarBatches{b};
        rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ctx.scenarioMaps{s}, ...
                                                               batch.rows, matRadCfg);
        oarCenter(batch.localIx, :) = oarCenter(batch.localIx, :) + ...
            scenarioWeight .* rows;
        if cacheOARRows
            blockBytes = ScenarioBatch.Cache.matRad_writeScenarioDoseCacheBlock(provider.cacheContext, s, ...
                                                                                oarCacheKind, b, batch.rows, rows);
            provider = ScenarioBatch.Resources.matRad_updateScenarioDoseDiskCacheUsage(provider, ...
                                                                                       blockBytes);
        end
        oarProgressReporter.update(sprintf( ...
                                           'scenario %d/%d, batch %d/%d', s, numScenarios, b, ...
                                           numel(oarBatches)));
    end

    source = [];
end

if ~isempty(ctx.targetRows)
    dijInterval.center(ctx.targetRows, :) = targetCenter;
end
if ~isempty(ctx.oarRows)
    dijInterval.center(ctx.oarRows, :) = oarCenter;
end
if strcmp(cfg.RadiusMode, 'std')
    dijInterval.radius = targetSecondMoment;
end
end

function result = matRad_computeIntervalScenarioBatchFirstPassScenario( ...
                                                                       provider, ctx, quantity, cfg, targetBatches, oarBatches, cacheTargetRows, ...
                                                                       cacheOARRows, scenarioIx, matRadCfg, targetProgressQueue, oarProgressQueue, ...
                                                                       numScenarios)
localProvider = ScenarioBatch.Worker.matRad_buildScenarioDoseWorkerProvider(provider);

targetCenter = sparse(numel(ctx.targetRows), ctx.numBixels);
oarCenter = sparse(numel(ctx.oarRows), ctx.numBixels);
targetSecondMoment = sparse(ctx.numBixels, ctx.numBixels);
targetCacheKind = 'target';
oarCacheKind = 'oar';

[localProvider, source] = ScenarioBatch.Source.matRad_createScenarioDoseRowSource( ...
                                                                                  localProvider, ctx, quantity, scenarioIx, matRadCfg);
scenarioWeight = ctx.scenarioWeights(scenarioIx);

for b = 1:numel(targetBatches)
    batch = targetBatches{b};
    rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ...
                                                           ctx.scenarioMaps{scenarioIx}, batch.rows, matRadCfg);
    targetCenter(batch.localIx, :) = targetCenter(batch.localIx, :) + ...
        scenarioWeight .* rows;
    if strcmp(cfg.RadiusMode, 'std')
        targetSecondMoment = targetSecondMoment + rows' * ...
            (scenarioWeight .* rows);
    end
    if cacheTargetRows
        blockBytes = ScenarioBatch.Cache.matRad_writeScenarioDoseCacheBlock( ...
                                                                            localProvider.cacheContext, scenarioIx, ...
                                                                            targetCacheKind, b, batch.rows, rows);
        localProvider = ScenarioBatch.Resources.matRad_updateScenarioDoseDiskCacheUsage( ...
                                                                                        localProvider, blockBytes);
    end
    if ~isempty(targetProgressQueue)
        send(targetProgressQueue, sprintf('scenario %d/%d, batch %d/%d', ...
                                          scenarioIx, numScenarios, b, numel(targetBatches)));
    end
end

for b = 1:numel(oarBatches)
    batch = oarBatches{b};
    rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ...
                                                           ctx.scenarioMaps{scenarioIx}, batch.rows, matRadCfg);
    oarCenter(batch.localIx, :) = oarCenter(batch.localIx, :) + ...
        scenarioWeight .* rows;
    if cacheOARRows
        blockBytes = ScenarioBatch.Cache.matRad_writeScenarioDoseCacheBlock( ...
                                                                            localProvider.cacheContext, scenarioIx, ...
                                                                            oarCacheKind, b, batch.rows, rows);
        localProvider = ScenarioBatch.Resources.matRad_updateScenarioDoseDiskCacheUsage( ...
                                                                                        localProvider, blockBytes);
    end
    if ~isempty(oarProgressQueue)
        send(oarProgressQueue, sprintf('scenario %d/%d, batch %d/%d', ...
                                       scenarioIx, numScenarios, b, numel(oarBatches)));
    end
end

source = [];

result = struct();
result.targetCenter = targetCenter;
result.oarCenter = oarCenter;
result.targetSecondMoment = targetSecondMoment;
result.resourceUsage = localProvider.resourceUsage;
end

function [dijInterval, provider] = matRad_accumulateIntervalScenarioBatchTargetExtremeRadius( ...
                                                                                             dijInterval, provider, ctx, quantity, ...
                                                                                             cfg, targetBatches, matRadCfg)
targetDeltaRows = sparse(numel(ctx.targetRows), ctx.numBixels);
numScenarios = numel(ctx.scenarioDijIx);

matRadCfg.dispInfo('matRad: Scenario-batch second pass for extreme target radius.\n');
stageName = 'scenario-batch interval target extreme radius';
resultBytesPerScenario = matRad_estimateIntervalTargetDeltaRowsBytes(ctx);
accumulatorBytes = 2 * resultBytesPerScenario;
[useParallel, parallelProvider, parallelPlan] = ScenarioBatch.Parallel.matRad_configureScenarioDoseParallel( ...
                                                                                                           provider, ctx, cfg, matRadCfg, ...
                                                                                                           stageName, [], ...
                                                                                                           resultBytesPerScenario, ...
                                                                                                           accumulatorBytes);
useParallel = DoseInterval.matRad_guardDoseIntervalTargetExtremeMemory(ctx, targetBatches, cfg, ...
                                                                       useParallel, matRadCfg, parallelPlan);
if useParallel
    logLevel = matRadCfg.logLevel;
    centerRowsAll = dijInterval.center(ctx.targetRows, :);
    progressWorkItems = numScenarios * numel(targetBatches);
    progressReporter = ScenarioBatch.Telemetry.matRad_createScenarioDoseProgressReporter( ...
                                                                                         matRadCfg, cfg, 'Target extreme radius', ...
                                                                                         progressWorkItems, true);
    progressCleanup = onCleanup(progressReporter.cleanup);
    progressQueue = progressReporter.queue;
    scenarioChunks = ScenarioBatch.Parallel.matRad_buildScenarioChunks( ...
                                                                       numScenarios, parallelPlan.chunkSize);
    for chunkIx = 1:numel(scenarioChunks)
        chunkScenarios = scenarioChunks{chunkIx};
        chunkResults = cell(numel(chunkScenarios), 1);
        parfor chunkLocalIx = 1:numel(chunkScenarios)
            scenarioIx = chunkScenarios(chunkLocalIx);
            workerCfg = MatRad_Config.instance();
            workerCfg.logLevel = logLevel;
            chunkResults{chunkLocalIx} = matRad_computeIntervalScenarioBatchTargetExtremeScenario( ...
                                                                                                  parallelProvider, ctx, quantity, cfg, targetBatches, ...
                                                                                                  centerRowsAll, ...
                                                                                                  scenarioIx, workerCfg, progressQueue, numScenarios);
        end

        for chunkLocalIx = 1:numel(chunkResults)
            result = chunkResults{chunkLocalIx};
            targetDeltaRows = max(targetDeltaRows, result.targetDeltaRows);
        end
        provider = ScenarioBatch.Resources.matRad_mergeScenarioDoseResourceUsage(provider, chunkResults);
        chunkResults = [];
    end
    if cfg.CollectTiming
        dijInterval.timing.parallelScenario.targetExtreme = true;
    end

    deltaSquared = full(sum(targetDeltaRows.^2, 1));
    centerRows = dijInterval.center(ctx.targetRows, :);
    dijInterval.radius = centerRows' * centerRows + ...
        spdiags(deltaSquared(:), 0, ctx.numBixels, ctx.numBixels);
    return
end

progressWorkItems = numScenarios * numel(targetBatches);
progressReporter = ScenarioBatch.Telemetry.matRad_createScenarioDoseProgressReporter( ...
                                                                                     matRadCfg, cfg, 'Target extreme radius', ...
                                                                                     progressWorkItems, false);
targetCacheKind = 'target';
for s = 1:numScenarios
    source = [];
    if strcmp(cfg.SecondPassStrategy, 'recompute')
        [provider, source] = ScenarioBatch.Source.matRad_createScenarioDoseRowSource( ...
                                                                                     provider, ctx, quantity, s, matRadCfg);
    end

    for b = 1:numel(targetBatches)
        batch = targetBatches{b};
        if strcmp(cfg.SecondPassStrategy, 'disk')
            rows = ScenarioBatch.Cache.matRad_readScenarioDoseCacheBlock(provider.cacheContext, s, ...
                                                                         targetCacheKind, b, batch.rows);
        else
            rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ctx.scenarioMaps{s}, ...
                                                                   batch.rows, matRadCfg);
        end
        centerRows = dijInterval.center(batch.rows, :);
        if strcmp(cfg.SecondPassStrategy, 'recompute')
            provider = ScenarioBatch.Resources.matRad_updateScenarioDoseMemoryUsage(provider, source, ...
                                                                                    rows, centerRows, targetDeltaRows);
        end
        targetDeltaRows(batch.localIx, :) = max(targetDeltaRows(batch.localIx, :), ...
                                                abs(rows - centerRows));
        progressReporter.update(sprintf('scenario %d/%d, batch %d/%d', ...
                                        s, numScenarios, b, numel(targetBatches)));
    end

    if strcmp(cfg.SecondPassStrategy, 'recompute')
        source = [];
    end
end

deltaSquared = full(sum(targetDeltaRows.^2, 1));
centerRows = dijInterval.center(ctx.targetRows, :);
dijInterval.radius = centerRows' * centerRows + ...
    spdiags(deltaSquared(:), 0, ctx.numBixels, ctx.numBixels);
end

function result = matRad_computeIntervalScenarioBatchTargetExtremeScenario( ...
                                                                           provider, ctx, quantity, cfg, targetBatches, centerRowsAll, scenarioIx, ...
                                                                           matRadCfg, progressQueue, numScenarios)
localProvider = ScenarioBatch.Worker.matRad_buildScenarioDoseWorkerProvider(provider);

targetDeltaRows = sparse(numel(ctx.targetRows), ctx.numBixels);
targetCacheKind = 'target';
source = [];
if strcmp(cfg.SecondPassStrategy, 'recompute')
    [localProvider, source] = ScenarioBatch.Source.matRad_createScenarioDoseRowSource( ...
                                                                                      localProvider, ctx, quantity, scenarioIx, matRadCfg);
end

for b = 1:numel(targetBatches)
    batch = targetBatches{b};
    if strcmp(cfg.SecondPassStrategy, 'disk')
        rows = ScenarioBatch.Cache.matRad_readScenarioDoseCacheBlock(localProvider.cacheContext, ...
                                                                     scenarioIx, targetCacheKind, b, batch.rows);
    else
        rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ...
                                                               ctx.scenarioMaps{scenarioIx}, batch.rows, matRadCfg);
    end
    centerRows = centerRowsAll(batch.localIx, :);
    if strcmp(cfg.SecondPassStrategy, 'recompute')
        localProvider = ScenarioBatch.Resources.matRad_updateScenarioDoseMemoryUsage(localProvider, ...
                                                                                     source, rows, centerRows, targetDeltaRows);
    end
    targetDeltaRows(batch.localIx, :) = max(targetDeltaRows(batch.localIx, :), ...
                                            abs(rows - centerRows));
    if ~isempty(progressQueue)
        send(progressQueue, sprintf('scenario %d/%d, batch %d/%d', ...
                                    scenarioIx, numScenarios, b, numel(targetBatches)));
    end
end

if strcmp(cfg.SecondPassStrategy, 'recompute')
    source = [];
end

result = struct();
result.targetDeltaRows = targetDeltaRows;
result.resourceUsage = localProvider.resourceUsage;
end

function [dijInterval, provider] = matRad_accumulateIntervalScenarioBatchOARRadiusFactors( ...
                                                                                          dijInterval, provider, ctx, quantity, ...
                                                                                          cfg, oarBatches, matRadCfg)
numScenarios = numel(ctx.scenarioDijIx);
oarRadiusOffset = 0;

matRadCfg.dispInfo('matRad: Scenario-batch second pass for INTERVAL3 OAR radius factors.\n');
progressReporter = ScenarioBatch.Telemetry.matRad_createScenarioDoseProgressReporter(matRadCfg, cfg, ...
                                                                                     'OAR radius factor', numel(oarBatches), false);
oarCacheKind = 'oar';
for b = 1:numel(oarBatches)
    batch = oarBatches{b};
    scenarioRows = cell(numScenarios, 1);
    progressReporter.update(sprintf('batch %d/%d', b, numel(oarBatches)));

    for s = 1:numScenarios
        if strcmp(cfg.SecondPassStrategy, 'disk')
            scenarioRows{s} = ScenarioBatch.Cache.matRad_readScenarioDoseCacheBlock(provider.cacheContext, ...
                                                                                    s, oarCacheKind, b, batch.rows);
        else
            [provider, source] = ScenarioBatch.Source.matRad_createScenarioDoseRowSource( ...
                                                                                         provider, ctx, quantity, s, matRadCfg);
            scenarioRows{s} = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ...
                                                                              ctx.scenarioMaps{s}, batch.rows, matRadCfg);
            provider = ScenarioBatch.Resources.matRad_updateScenarioDoseMemoryUsage(provider, source, ...
                                                                                    scenarioRows{s});
            source = [];
        end
    end
    if strcmp(cfg.SecondPassStrategy, 'recompute')
        provider = ScenarioBatch.Resources.matRad_updateScenarioDoseMemoryUsage(provider, scenarioRows);
    end

    for localIx = 1:numel(batch.rows)
        centeredMatrixRows = cell(numScenarios, 1);
        centerRow = dijInterval.center(batch.rows(localIx), :);
        for s = 1:numScenarios
            centeredMatrixRows{s} = scenarioRows{s}(localIx, :) - centerRow;
        end
        centeredScenarioMatrix = vertcat(centeredMatrixRows{:});
        if strcmp(cfg.SecondPassStrategy, 'recompute')
            provider = ScenarioBatch.Resources.matRad_updateScenarioDoseMemoryUsage(provider, ...
                                                                                    scenarioRows, centerRow, centeredScenarioMatrix);
        end
        scenarioWeights = ctx.scenarioWeights;
        [radiusFactor, rank] = DoseInterval.matRad_buildDoseIntervalRadiusFactor( ...
                                                                                 centeredScenarioMatrix, scenarioWeights, cfg, ctx.numBixels);
        dijInterval.OARRadiusRank(oarRadiusOffset + localIx) = rank;
        dijInterval.OARRadiusFactor{oarRadiusOffset + localIx} = radiusFactor;
    end

    oarRadiusOffset = oarRadiusOffset + numel(batch.rows);
end
end

function bytes = matRad_estimateIntervalFirstPassResultBytes(ctx, cfg)
sparseFillFactor = 0.05;
targetCenterBytes = ScenarioBatch.Resources.matRad_estimateSparseMatrixBytes( ...
                                                                             numel(ctx.targetRows), ctx.numBixels, sparseFillFactor);
oarCenterBytes = ScenarioBatch.Resources.matRad_estimateSparseMatrixBytes( ...
                                                                          numel(ctx.oarRows), ctx.numBixels, sparseFillFactor);
targetSecondMomentBytes = 0;
if strcmp(cfg.RadiusMode, 'std')
    targetSecondMomentBytes = double(ctx.numBixels)^2 * 8;
end
bytes = targetCenterBytes + oarCenterBytes + targetSecondMomentBytes;
end

function bytes = matRad_estimateIntervalTargetDeltaRowsBytes(ctx)
sparseFillFactor = 0.05;
bytes = ScenarioBatch.Resources.matRad_estimateSparseMatrixBytes( ...
                                                                 numel(ctx.targetRows), ctx.numBixels, sparseFillFactor);
end
