function [pln_prob2,dij_prob2Context] = matRad_calcDoseProb2(ct,cst,stf,pln,dij,cfg)
% matRad_calcDoseProb2 calculates scenario-free probabilistic influence data
%
% call
%   [pln_prob2,dij_prob2Context] = matRad_calcDoseProb2(ct,cst,stf,pln,dij)
%   [pln_prob2,dij_prob2Context] = matRad_calcDoseProb2(ct,cst,stf,pln,dij,cfg)
%
% input
%   ct:     matRad ct struct
%   cst:    matRad cst cell array
%   stf:    matRad steering information struct (reserved)
%   pln:    matRad pln struct with robust scenario model; optional
%           pln.propOpt.scen4D selects CT scenarios used for the
%           calculation (default: 1; 'all' or positive integer CT scenario
%           ids)
%   dij:    matRad dose influence struct for the robust scenarios
%   cfg:    optional configuration struct with fields compatible with
%           matRad_calcDoseInterval2 for Quantity, QuantityField, refScen,
%           targetStructSel, OARStructSel, MemoryLimitMB, BatchSize, and
%           ProgressLevel
%
% output
%   pln_prob2:         plan struct; pln_prob2.propOpt.dij_prob2 contains
%                      selected-VOI expected dose influence and VOI
%                      variance matrices
%   dij_prob2Context:  lightweight dij context passed to
%                      matRad_fluenceOptimization
%
% References
%   Cristoforetti et al., Med Phys 2025, scenario-free robust optimization.
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

if nargin < 6
    cfg = struct();
end

matRad_cfg = MatRad_Config.instance();
timer = tic;

ctx = matRad_resolveScenarioDoseInputs(ct,cst,pln,dij,cfg, ...
    'PROB2',matRad_cfg);
cfg = ctx.cfg;
quantity = validateProb2Quantity(ctx.quantity,matRad_cfg);

matRad_cfg.dispInfo(['matRad: Calculating PROB2 scenario-free data for ', ...
    'quantity ''%s'' using dij.%s and %d scenarios.\n'], ...
    quantity.name,quantity.field,numel(ctx.scenarioDijIx));
matRad_cfg.dispInfo(['matRad: PROB2 reference CT scenario %d, %d voxels, ', ...
    '%d bixels.\n'],cfg.refScen,ctx.numVoxels,ctx.numBixels);

dij_prob2 = initializeProb2Struct(ctx,quantity);
expectedRows = resolveProb2ExpectedRows(ctx);
dij_prob2.expected = accumulateExpectedInfluence(quantity,ctx,cfg, ...
    expectedRows,matRad_cfg);

voiRows = resolveProb2VoiRows(ctx);
dij_prob2.voiSubIx = voiRows;
dij_prob2.Omega = accumulateVoiOmega(quantity,ctx,cfg,voiRows,matRad_cfg);

pln_prob2 = pln;
if ~isfield(pln_prob2,'propOpt') || ~isstruct(pln_prob2.propOpt)
    pln_prob2.propOpt = struct();
end
pln_prob2.propOpt.dij_prob2 = dij_prob2;
dij_prob2Context = buildProb2DijContext(ct,cst,stf,pln,dij, ...
    dij_prob2,quantity,cfg);
pln_prob2.multScen = dij_prob2Context.scenarioModel;

matRad_cfg.dispInfo('matRad: Finished PROB2 calculation in %.2f s.\n', ...
    toc(timer));

end

function quantity = validateProb2Quantity(quantity,matRad_cfg)
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

function dij_prob2 = initializeProb2Struct(ctx,quantity)
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

function expected = accumulateExpectedInfluence(quantity,ctx,cfg,rowsAll, ...
    matRad_cfg)
numVoxels = ctx.numVoxels;
numBixels = ctx.numBixels;
batchSize = resolveProb2BatchSize(numel(rowsAll), ...
    numel(ctx.scenarioDijIx),numBixels,cfg);
batches = makeProb2Batches(rowsAll,batchSize);
expected = sparse(numVoxels,numBixels);

matRad_cfg.dispInfo(['matRad: Accumulating PROB2 expected influence for ', ...
    '%d selected voxels in %d batches of up to %d voxels.\n'], ...
    numel(rowsAll),numel(batches),batchSize);

for b = 1:numel(batches)
    rows = batches{b};
    logProb2BatchProgress(matRad_cfg,cfg,'Expected influence',b,numel(batches));
    options.StoreScenarioRows = false;
    [stats,~] = computeProb2ScenarioDoseBatchStats(quantity,ctx,rows, ...
        matRad_cfg,options);
    expected(rows,:) = stats.meanRows;
end
end

function [stats,timing] = computeProb2ScenarioDoseBatchStats(quantity,ctx, ...
    rows,matRad_cfg,options)
[stats,timing] = matRad_computeScenarioDoseBatchStats(quantity, ...
    ctx.scenarioDijIx,ctx.scenarioWeights,ctx.scenarioMaps,rows, ...
    ctx.numBixels,matRad_cfg,options);
end

function Omega = accumulateVoiOmega(quantity,ctx,cfg,voiRows,matRad_cfg)
Omega = cell(size(voiRows));
numBixels = ctx.numBixels;

for i = 1:numel(voiRows)
    rows = voiRows{i};
    if isempty(rows)
        Omega{i} = [];
        continue;
    end

    matRad_cfg.dispInfo(['matRad: Accumulating PROB2 Omega for VOI %d ', ...
        '(%d voxels).\n'],i,numel(rows));

    batchSize = resolveProb2BatchSize(numel(rows), ...
        numel(ctx.scenarioDijIx),numBixels,cfg);
    batches = makeProb2Batches(rows,batchSize);
    omega = sparse(numBixels,numBixels);
    options.ComputeCenteredRows = true;

    for b = 1:numel(batches)
        batchRows = batches{b};
        logProb2BatchProgress(matRad_cfg,cfg,sprintf('Omega VOI %d',i), ...
            b,numel(batches));

        [stats,~] = computeProb2ScenarioDoseBatchStats(quantity,ctx, ...
            batchRows,matRad_cfg,options);

        omega = omega + accumulateCenteredOmegaBatch(stats.centeredRows, ...
            stats.scenarioWeights,numBixels);
    end

    Omega{i} = sparse(0.5.*(omega + omega'));
end
end

function omega = accumulateCenteredOmegaBatch(centeredRows,scenarioWeights, ...
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

function rows = resolveProb2ExpectedRows(ctx)
rows = unique([ctx.targetRows(:); ctx.oarRows(:)],'stable');
end

function voiRows = resolveProb2VoiRows(ctx)
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

function dij_prob2Context = buildProb2DijContext(ct,cst,stf,pln,dij, ...
    dij_prob2,quantity,cfg)
numBixels = size(dij_prob2.expected,2);
dij_prob2Context = matRad_buildNominalDijContext( ...
    ct,cst,stf,pln,cfg,numBixels,dij);
if strcmpi(quantity.name,'RBExD') && isfield(dij_prob2Context,'RBE')
    dij_prob2Context.RBE = 1;
end
dij_prob2Context.probabilisticQuantity = quantity.name;
dij_prob2Context.probabilisticQuantityField = quantity.field;
end

function batchSize = resolveProb2BatchSize(numRows,numScenarios,numBixels,cfg)
if numRows == 0
    batchSize = 0;
    return;
end

if ~isempty(cfg.BatchSize)
    batchSize = min(numRows,cfg.BatchSize);
    return;
end

memoryLimitMB = cfg.MemoryLimitMB;
if isempty(memoryLimitMB)
    memoryLimitMB = 256;
end
bytesPerRow = max(1,numScenarios*numBixels*8);
batchSize = max(1,floor(memoryLimitMB*1e6/bytesPerRow));
batchSize = min(numRows,batchSize);
end

function batches = makeProb2Batches(rows,batchSize)
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

function logProb2BatchProgress(matRad_cfg,cfg,stageName,batchIx,numBatches)
if numBatches == 0
    return;
end

progressStep = max(1,ceil(numBatches/10));
if isfield(cfg,'ProgressLevel') && ...
        (strcmp(cfg.ProgressLevel,'detailed') || batchIx == 1 || ...
        batchIx == numBatches || mod(batchIx,progressStep) == 0)
    matRad_cfg.dispInfo('matRad: %s progress %d/%d (%.0f%%).\n', ...
        stageName,batchIx,numBatches,100*batchIx/numBatches);
end
end
