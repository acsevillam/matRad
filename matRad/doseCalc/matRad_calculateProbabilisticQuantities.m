function [plnProb, dijProbContext] = matRad_calculateProbabilisticQuantities(ct, cst, stf, pln, varargin)
% matRad_calculateProbabilisticQuantities calculates probabilistic data scenario by scenario
%
% call
%   [plnProb,dijProbContext] = matRad_calculateProbabilisticQuantities(ct,cst,stf,pln)
%   [plnProb,dijProbContext] = matRad_calculateProbabilisticQuantities(ct,cst,stf,pln,dij)
%   [plnProb,dijProbContext] = matRad_calculateProbabilisticQuantities(ct,cst,stf,pln,cfg)
%   [plnProb,dijProbContext] = matRad_calculateProbabilisticQuantities(ct,cst,stf,pln,dij,cfg)
%
% input
%   ct:     matRad ct struct
%   cst:    matRad cst cell array
%   stf:    matRad steering information struct
%   pln:    matRad pln struct with robust scenario model
%   dij:    optional precomputed robust dose influence struct. If provided,
%           precomputed scenario matrices from dij are used instead of
%           recalculating scenario dose influence data.
%   cfg:    optional configuration struct with fields:
%           UseParallel: use safe available scenario parallelism for the
%               first pass and omega pass when the Parallel Computing
%               Toolbox and enough workers/memory are available. If omega
%               parallel aggregation exceeds MemoryLimitMB but serial omega
%               fits, only the omega pass falls back to serial execution.
%               matRad may create, reduce, or increase the active parallel
%               pool.
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
%   plnProb:        plan struct containing selected-VOI propOpt.dij_prob.
%                     The propOpt.dij_prob.precomputeSize field summarizes
%                     compact result bytes and peak auxiliary bytes used
%                     during precomputation.
%   dijProbContext: lightweight single-scenario dij context for
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

matRadCfg = MatRad_Config.instance();
[~, cfg] = ScenarioBatch.Config.matRad_parseScenarioDosePrecomputeArguments(matRadCfg, varargin{:});
timer = tic;
cacheContext = [];

try
    cfg = DoseProb.matRad_normalizeDoseProbConfig(cfg, matRadCfg);
    scenarioInfo = ScenarioBatch.Input.matRad_selectScenarioDoseRealizations(ct, pln, cfg, matRadCfg);
    [provider, dijForResolve] = ScenarioBatch.Source.matRad_createScenarioDoseRowProvider( ...
                                                                                          ct, cst, stf, pln, cfg, ...
                                                                                          scenarioInfo, matRadCfg, 'PROB');
    provider.resourceUsage = ScenarioBatch.Resources.matRad_buildScenarioDoseResourceUsage();

    ctx = ScenarioBatch.Input.matRad_buildScenarioDoseContext(ct, cst, pln, dijForResolve, cfg, ...
                                                              'PROB', matRadCfg, scenarioInfo);
    ctx.scenarioIds = scenarioInfo.scenarioIds(:);
    if isfield(dijForResolve, 'doseGrid')
        ctx.doseGrid = dijForResolve.doseGrid;
    end
    if isfield(dijForResolve, 'ctGrid')
        ctx.ctGrid = dijForResolve.ctGrid;
    end

    cfg = ctx.cfg;
    cfg.SecondPassStrategy = provider.secondPassStrategy;
    cfg.KeepCache = provider.keepCache;
    cfg.CacheRoot = provider.cacheRoot;
    quantity = ctx.quantity;

    matRadCfg.dispInfo(['matRad: Calculating probabilistic scenario-batch data for ', ...
                        'quantity ''%s'' and %d scenarios.\n'], ...
                       quantity.name, numel(ctx.scenarioDijIx));
    matRadCfg.dispInfo(['matRad: Probabilistic scenario-batch reference CT scenario ', ...
                        '%d, %d voxels, %d bixels.\n'], ...
                       cfg.refScen, ctx.numVoxels, ctx.numBixels);

    [~, expectedBatches, voiRows, voiBatches] = DoseProb.matRad_buildDoseProbBatches(ctx, cfg);

    needsSecondPass = any(~cellfun(@isempty, voiRows(:)));
    needsDiskCache = strcmp(cfg.SecondPassStrategy, 'disk') && needsSecondPass;
    if needsDiskCache
        signatureExtras = struct();
        signatureExtras.voiRows = voiRows;
        cacheContext = ScenarioBatch.Cache.matRad_createScenarioDoseCache(cfg, ctx, quantity, ...
                                                                          stf, matRadCfg, 'PROB', signatureExtras);
        provider.cacheContext = cacheContext;
    end

    firstPass = struct();
    firstPass.expectedBatches = expectedBatches;
    firstPass.voiBatches = voiBatches;
    firstPass.cacheVoiRows = needsDiskCache;

    dijProb = DoseProb.matRad_buildDoseProbDij(ctx, quantity);
    dijProb.voiSubIx = voiRows;
    dijProb.precomputeMode = 'scenario-batch';
    dijProb.secondPassStrategy = cfg.SecondPassStrategy;

    [dijProb, provider] = matRad_accumulateProbScenarioBatchFirstPass( ...
                                                                      dijProb, provider, ctx, quantity, cfg, firstPass, matRadCfg);

    if needsSecondPass
        [dijProb.Omega, provider, omegaParallel] = matRad_accumulateProbScenarioBatchOmega( ...
                                                                                           dijProb, provider, ctx, quantity, ...
                                                                                           cfg, voiRows, voiBatches, matRadCfg);
        if cfg.CollectTiming
            dijProb.timing.parallelScenario.omega = omegaParallel;
        end
    else
        dijProb.Omega = cell(size(voiRows));
    end

    if cfg.CollectTiming
        dijProb.timing.totalSeconds = toc(timer);
    end

    if needsDiskCache && cfg.KeepCache
        dijProb.cacheDir = cacheContext.runDir;
    end
    dijProb.precomputeSize = ScenarioBatch.Resources.matRad_buildScenarioDosePrecomputeSize( ...
                                                                                            dijProb, provider, cfg);

    plnProb = pln;
    if ~isfield(plnProb, 'propOpt') || ~isstruct(plnProb.propOpt)
        plnProb.propOpt = struct();
    end
    plnProb.propOpt.dij_prob = dijProb;
    dijProbContext = DoseProb.matRad_buildDoseProbDijContext(ct, cst, stf, pln, ...
                                                             dijProb, quantity, cfg);
    plnProb.multScen = dijProbContext.scenarioModel;

    ScenarioBatch.Cache.matRad_cleanupScenarioDoseCache(cacheContext, cfg.KeepCache);
    matRadCfg.dispInfo('matRad: Finished probabilistic scenario-batch calculation in %.2f s.\n', ...
                       toc(timer));
catch ME
    if ~isempty(cacheContext)
        ScenarioBatch.Cache.matRad_cleanupScenarioDoseCache(cacheContext, false);
    end
    rethrow(ME);
end

end

function [dijProb, provider] = matRad_accumulateProbScenarioBatchFirstPass( ...
                                                                           dijProb, provider, ctx, quantity, cfg, firstPass, matRadCfg)
numScenarios = numel(ctx.scenarioDijIx);
expectedBatches = firstPass.expectedBatches;
voiBatches = firstPass.voiBatches;
cacheVoiRows = firstPass.cacheVoiRows;
numExpectedRows = sum(cellfun(@(batch) numel(batch.rows), expectedBatches));

matRadCfg.dispInfo(['matRad: Scenario-batch probabilistic first pass over %d ', ...
                    'scenario(s) and %d selected dose voxels.\n'], ...
                   numScenarios, numExpectedRows);

stageName = 'scenario-batch probabilistic first pass';
[useParallel, parallelProvider] = ScenarioBatch.Parallel.matRad_configureScenarioDoseParallel( ...
                                                                                              provider, ctx, cfg, matRadCfg, stageName, []);
if useParallel
    expectedRows = sparse(numExpectedRows, ctx.numBixels);
    logLevel = matRadCfg.logLevel;
    progressWorkItems = numScenarios * numel(expectedBatches);
    progressReporter = ScenarioBatch.Telemetry.matRad_createScenarioDoseProgressReporter( ...
                                                                                         matRadCfg, cfg, 'Expected influence', ...
                                                                                         progressWorkItems, true);
    progressCleanup = onCleanup(progressReporter.cleanup);
    progressQueue = progressReporter.queue;
    scenarioResults = cell(numScenarios, 1);
    parfor s = 1:numScenarios
        workerCfg = MatRad_Config.instance();
        workerCfg.logLevel = logLevel;
        scenarioResults{s} = matRad_computeProbScenarioBatchFirstPassScenario( ...
                                                                              parallelProvider, ctx, quantity, cfg, expectedBatches, voiBatches, ...
                                                                              cacheVoiRows, s, numExpectedRows, workerCfg, progressQueue, ...
                                                                              numScenarios);
    end

    for s = 1:numScenarios
        result = scenarioResults{s};
        expectedRows = expectedRows + result.expectedRows;
    end
    provider = ScenarioBatch.Resources.matRad_mergeScenarioDoseResourceUsage(provider, scenarioResults);
    if cfg.CollectTiming
        dijProb.timing.parallelScenario.firstPass = true;
    end

    for b = 1:numel(expectedBatches)
        batch = expectedBatches{b};
        dijProb.expected(batch.rows, :) = expectedRows(batch.localIx, :);
    end
    return
end

progressWorkItems = numScenarios * numel(expectedBatches);
progressReporter = ScenarioBatch.Telemetry.matRad_createScenarioDoseProgressReporter( ...
                                                                                     matRadCfg, cfg, 'Expected influence', ...
                                                                                     progressWorkItems, false);
for s = 1:numScenarios
    [provider, source] = ScenarioBatch.Source.matRad_createScenarioDoseRowSource(provider, ctx, ...
                                                                                 quantity, s, matRadCfg);
    scenarioWeight = ctx.scenarioWeights(s);

    for b = 1:numel(expectedBatches)
        batch = expectedBatches{b};
        rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ctx.scenarioMaps{s}, ...
                                                               batch.rows, matRadCfg);
        dijProb.expected(batch.rows, :) = dijProb.expected(batch.rows, :) + ...
            scenarioWeight .* rows;
        progressReporter.update(sprintf('scenario %d/%d, batch %d/%d', ...
                                        s, numScenarios, b, numel(expectedBatches)));
    end

    if cacheVoiRows
        for voiIx = 1:numel(voiBatches)
            batches = voiBatches{voiIx};
            cacheKind = sprintf('voi_%04d', voiIx);
            for b = 1:numel(batches)
                batch = batches{b};
                rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ...
                                                                       ctx.scenarioMaps{s}, batch.rows, matRadCfg);
                blockBytes = ScenarioBatch.Cache.matRad_writeScenarioDoseCacheBlock( ...
                                                                                    provider.cacheContext, s, cacheKind, b, batch.rows, rows);
                provider = ScenarioBatch.Resources.matRad_updateScenarioDoseDiskCacheUsage(provider, ...
                                                                                           blockBytes);
            end
        end
    end

    source = [];
end
end

function result = matRad_computeProbScenarioBatchFirstPassScenario( ...
                                                                   provider, ctx, quantity, cfg, expectedBatches, voiBatches, cacheVoiRows, ...
                                                                   scenarioIx, numExpectedRows, matRadCfg, progressQueue, numScenarios)
localProvider = ScenarioBatch.Worker.matRad_buildScenarioDoseWorkerProvider(provider);

expectedRows = sparse(numExpectedRows, ctx.numBixels);
[localProvider, source] = ScenarioBatch.Source.matRad_createScenarioDoseRowSource( ...
                                                                                  localProvider, ctx, quantity, scenarioIx, matRadCfg);
scenarioWeight = ctx.scenarioWeights(scenarioIx);

for b = 1:numel(expectedBatches)
    batch = expectedBatches{b};
    rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ...
                                                           ctx.scenarioMaps{scenarioIx}, batch.rows, matRadCfg);
    expectedRows(batch.localIx, :) = expectedRows(batch.localIx, :) + ...
        scenarioWeight .* rows;
    if ~isempty(progressQueue)
        send(progressQueue, sprintf('scenario %d/%d, batch %d/%d', ...
                                    scenarioIx, numScenarios, b, numel(expectedBatches)));
    end
end

if cacheVoiRows
    for voiIx = 1:numel(voiBatches)
        batches = voiBatches{voiIx};
        cacheKind = sprintf('voi_%04d', voiIx);
        for b = 1:numel(batches)
            batch = batches{b};
            rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ...
                                                                   ctx.scenarioMaps{scenarioIx}, batch.rows, matRadCfg);
            blockBytes = ScenarioBatch.Cache.matRad_writeScenarioDoseCacheBlock( ...
                                                                                localProvider.cacheContext, scenarioIx, ...
                                                                                cacheKind, b, batch.rows, rows);
            localProvider = ScenarioBatch.Resources.matRad_updateScenarioDoseDiskCacheUsage( ...
                                                                                            localProvider, blockBytes);
        end
    end
end

source = [];

result = struct();
result.expectedRows = expectedRows;
result.resourceUsage = localProvider.resourceUsage;
end

function [omegaCells, provider, omegaParallel] = matRad_accumulateProbScenarioBatchOmega( ...
                                                                                         dijProb, provider, ctx, quantity, cfg, ...
                                                                                         voiRows, voiBatches, matRadCfg)
omegaCells = cell(size(voiRows));
numScenarios = numel(ctx.scenarioDijIx);
numBixels = ctx.numBixels;
omegaParallel = false;

stageName = 'scenario-batch probabilistic omega';
[useParallel, parallelProvider] = ScenarioBatch.Parallel.matRad_configureScenarioDoseParallel( ...
                                                                                              provider, ctx, cfg, matRadCfg, stageName, []);
useParallel = DoseProb.matRad_guardDoseProbOmegaMemory(ctx, voiBatches, cfg, useParallel, matRadCfg);
if useParallel
    logLevel = matRadCfg.logLevel;
    numOmegaBatches = sum(cellfun(@numel, voiBatches));
    progressReporter = ScenarioBatch.Telemetry.matRad_createScenarioDoseProgressReporter(matRadCfg, ...
                                                                                         cfg, 'omega', numScenarios * numOmegaBatches, true);
    progressCleanup = onCleanup(progressReporter.cleanup);
    progressQueue = progressReporter.queue;
    scenarioResults = cell(numScenarios, 1);
    parfor s = 1:numScenarios
        workerCfg = MatRad_Config.instance();
        workerCfg.logLevel = logLevel;
        scenarioResults{s} = matRad_computeProbScenarioBatchOmegaScenario( ...
                                                                          parallelProvider, dijProb, ctx, quantity, cfg, voiBatches, s, ...
                                                                          workerCfg, progressQueue, numScenarios);
    end

    for voiIx = 1:numel(voiRows)
        if isempty(voiRows{voiIx})
            omegaCells{voiIx} = [];
        else
            omegaCells{voiIx} = sparse(numBixels, numBixels);
        end
    end
    for s = 1:numScenarios
        result = scenarioResults{s};
        for voiIx = 1:numel(voiRows)
            if ~isempty(result.omegaCells{voiIx})
                omegaCells{voiIx} = omegaCells{voiIx} + result.omegaCells{voiIx};
            end
        end
    end
    provider = ScenarioBatch.Resources.matRad_mergeScenarioDoseResourceUsage(provider, scenarioResults);
    for voiIx = 1:numel(omegaCells)
        if ~isempty(omegaCells{voiIx})
            omegaCells{voiIx} = sparse(0.5 .* (omegaCells{voiIx} + omegaCells{voiIx}'));
        end
    end
    omegaParallel = true;
    return
end

numOmegaBatches = sum(cellfun(@numel, voiBatches));
progressReporter = ScenarioBatch.Telemetry.matRad_createScenarioDoseProgressReporter(matRadCfg, ...
                                                                                     cfg, 'omega', numScenarios * numOmegaBatches, false);
for voiIx = 1:numel(voiRows)
    if isempty(voiRows{voiIx})
        omegaCells{voiIx} = [];
        continue
    end

    matRadCfg.dispInfo(['matRad: Scenario-batch probabilistic second pass omega for ', ...
                        'VOI %d (%d voxels).\n'], voiIx, numel(voiRows{voiIx}));
    omegaMatrix = sparse(numBixels, numBixels);
    batches = voiBatches{voiIx};
    cacheKind = sprintf('voi_%04d', voiIx);

    if strcmp(cfg.SecondPassStrategy, 'disk')
        for b = 1:numel(batches)
            batch = batches{b};
            centeredRows = cell(numScenarios, 1);
            for s = 1:numScenarios
                rows = ScenarioBatch.Cache.matRad_readScenarioDoseCacheBlock(provider.cacheContext, ...
                                                                             s, cacheKind, b, batch.rows);
                centeredRows{s} = rows - dijProb.expected(batch.rows, :);
                progressReporter.update(sprintf( ...
                                                'scenario %d/%d, VOI %d, batch %d/%d', ...
                                                s, numScenarios, voiIx, b, numel(batches)));
            end
            omegaBlock = DoseProb.matRad_buildDoseProbOmegaBlock( ...
                                                                 centeredRows, ctx.scenarioWeights, numBixels);
            omegaMatrix = omegaMatrix + omegaBlock;
        end
    else
        for s = 1:numScenarios
            [provider, source] = ScenarioBatch.Source.matRad_createScenarioDoseRowSource(provider, ...
                                                                                         ctx, quantity, s, matRadCfg);
            for b = 1:numel(batches)
                batch = batches{b};
                rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ...
                                                                       ctx.scenarioMaps{s}, batch.rows, matRadCfg);
                centeredRows = {rows - dijProb.expected(batch.rows, :)};
                provider = ScenarioBatch.Resources.matRad_updateScenarioDoseMemoryUsage(provider, ...
                                                                                        source, rows, centeredRows);
                scenarioWeight = ctx.scenarioWeights(s);
                omegaBlock = DoseProb.matRad_buildDoseProbOmegaBlock( ...
                                                                     centeredRows, scenarioWeight, numBixels);
                omegaMatrix = omegaMatrix + omegaBlock;
                progressReporter.update(sprintf( ...
                                                'scenario %d/%d, VOI %d, batch %d/%d', ...
                                                s, numScenarios, voiIx, b, numel(batches)));
            end
            source = [];
        end
    end

    omegaCells{voiIx} = sparse(0.5 .* (omegaMatrix + omegaMatrix'));
end
end

function result = matRad_computeProbScenarioBatchOmegaScenario( ...
                                                               provider, dijProb, ctx, quantity, cfg, voiBatches, scenarioIx, matRadCfg, ...
                                                               progressQueue, numScenarios)
localProvider = ScenarioBatch.Worker.matRad_buildScenarioDoseWorkerProvider(provider);

omegaCells = cell(size(dijProb.voiSubIx));
numBixels = ctx.numBixels;
source = [];
if strcmp(cfg.SecondPassStrategy, 'recompute')
    [localProvider, source] = ScenarioBatch.Source.matRad_createScenarioDoseRowSource( ...
                                                                                      localProvider, ctx, quantity, scenarioIx, matRadCfg);
end

for voiIx = 1:numel(voiBatches)
    batches = voiBatches{voiIx};
    if isempty(batches)
        omegaCells{voiIx} = [];
        continue
    end

    omegaMatrix = sparse(numBixels, numBixels);
    cacheKind = sprintf('voi_%04d', voiIx);
    for b = 1:numel(batches)
        batch = batches{b};
        if strcmp(cfg.SecondPassStrategy, 'disk')
            rows = ScenarioBatch.Cache.matRad_readScenarioDoseCacheBlock(localProvider.cacheContext, ...
                                                                         scenarioIx, cacheKind, b, batch.rows);
        else
            rows = ScenarioBatch.Source.matRad_getScenarioDoseRows(source, ...
                                                                   ctx.scenarioMaps{scenarioIx}, batch.rows, matRadCfg);
        end
        centeredRows = {rows - dijProb.expected(batch.rows, :)};
        if strcmp(cfg.SecondPassStrategy, 'recompute')
            localProvider = ScenarioBatch.Resources.matRad_updateScenarioDoseMemoryUsage( ...
                                                                                         localProvider, source, rows, centeredRows);
        end
        scenarioWeight = ctx.scenarioWeights(scenarioIx);
        omegaBlock = DoseProb.matRad_buildDoseProbOmegaBlock( ...
                                                             centeredRows, scenarioWeight, numBixels);
        omegaMatrix = omegaMatrix + omegaBlock;
        if ~isempty(progressQueue)
            send(progressQueue, sprintf( ...
                                        'scenario %d/%d, VOI %d, batch %d/%d', ...
                                        scenarioIx, numScenarios, voiIx, b, numel(batches)));
        end
    end
    omegaCells{voiIx} = omegaMatrix;
end

if strcmp(cfg.SecondPassStrategy, 'recompute')
    source = [];
end

result = struct();
result.omegaCells = omegaCells;
result.resourceUsage = localProvider.resourceUsage;
end
