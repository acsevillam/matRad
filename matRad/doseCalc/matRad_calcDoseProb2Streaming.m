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
%           SecondPassStrategy: 'disk' (default) or 'recompute'
%           CacheRoot: root folder for disk blocks
%           KeepCache: keep the hash cache folder after the run (default false)
%
% output
%   pln_prob2:        plan struct containing propOpt.dij_prob2. The
%                     propOpt.dij_prob2.streamingSize field summarizes
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

    rowsAll = (1:ctx.numVoxels)';
    allBatches = matRad_makeScenarioDoseRowBatches(rowsAll, ...
        matRad_resolveScenarioDoseBatchSize(numel(rowsAll), ...
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
    firstPassWork.primaryBatches = allBatches;
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
        [dij_prob2.Omega,provider] = accumulateProb2StreamingOmega( ...
            dij_prob2,provider,ctx,quantity,cfg,voiRows,voiBatches, ...
            matRad_cfg);
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

matRad_cfg.dispInfo(['matRad: Streaming PROB2 first pass over %d ', ...
    'scenario(s) and %d dose voxels.\n'],numScenarios,ctx.numVoxels);

for s = 1:numScenarios
    matRad_cfg.dispInfo('matRad: Streaming PROB2 first pass scenario %d/%d.\n', ...
        s,numScenarios);
    [provider,source] = matRad_beginScenarioDoseRowsProvider(provider,ctx, ...
        quantity,s,matRad_cfg);
    scenarioWeight = ctx.scenarioWeights(s);

    for b = 1:numel(expectedBatches)
        batch = expectedBatches{b};
        matRad_logScenarioDoseBatchProgress(matRad_cfg,cfg,'Expected influence', ...
            b,numel(expectedBatches));
        rows = matRad_getScenarioDoseProviderRows(source,ctx.scenarioMaps{s}, ...
            batch.rows,matRad_cfg);
        dij_prob2.expected(batch.rows,:) = dij_prob2.expected(batch.rows,:) + ...
            scenarioWeight .* rows;
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

function [Omega,provider] = accumulateProb2StreamingOmega( ...
    dij_prob2,provider,ctx,quantity,cfg,voiRows,voiBatches,matRad_cfg)
Omega = cell(size(voiRows));
numScenarios = numel(ctx.scenarioDijIx);
numBixels = ctx.numBixels;

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
            end
            [provider,source] = matRad_endScenarioDoseRowsProvider(provider,source);
        end
    end

    Omega{voiIx} = sparse(0.5.*(omega + omega'));
end
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
