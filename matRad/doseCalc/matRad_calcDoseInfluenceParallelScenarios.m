function [dij,useParallel] = matRad_calcDoseInfluenceParallelScenarios(ct,cst,stf,pln,engine)
% matRad_calcDoseInfluenceParallelScenarios computes robust pencil-beam dijs in parallel
%
% call
%   [dij,useParallel] = matRad_calcDoseInfluenceParallelScenarios(ct,cst,stf,pln,engine)
%
% input
%   ct:         ct cube
%   cst:        matRad cst cell array
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   engine:     configured matRad dose engine
%
% output
%   dij:            assembled robust dij struct, or [] if serial fallback
%                   should be used
%   useParallel:    true if a parallel scenario calculation was executed
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
dij = [];
useParallel = false;

if nargin < 5 || isempty(engine)
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
end

if ~isa(engine,'DoseEngines.matRad_PencilBeamEngineAbstract')
    matRad_cfg.dispWarning(['UseParallel was requested for multi-scenario ', ...
        'dij calculation, but dose engine "%s" is not an analytical ', ...
        'pencil-beam engine. Falling back to serial dose calculation.\n'], ...
        engine.name);
    return;
end

scenarioModel = resolveScenarioModel(ct,pln,engine,matRad_cfg);
if isempty(scenarioModel)
    return;
end

scenarioIds = scenarioModel.scenarioIds();
numScenarios = numel(scenarioIds);
if numScenarios < 2
    matRad_cfg.dispWarning(['UseParallel was requested for multi-scenario ', ...
        'dij calculation, but fewer than two scenarios are active. ', ...
        'Falling back to serial dose calculation.\n']);
    return;
end

plnParallel = matRad_makeWorkerSafePlan(pln);
plnParallel.multScen = scenarioModel;

workerMemoryBytes = estimateScenarioDoseWorkerMemoryBytes(ct,cst,stf, ...
    plnParallel,engine);
[poolReady,~,~] = matRad_configureSafeDoseParallelPool( ...
    workerMemoryBytes,numScenarios,matRad_cfg, ...
    'multi-scenario dij calculation', ...
    'fallbackDescription','serial dose calculation');
if ~poolReady
    return;
end

scenarioDijs = cell(numScenarios,1);

parfor scenarioIx = 1:numScenarios
    scenarioId = scenarioIds(scenarioIx);
    plnScenario = matRad_prepareSerialScenarioPlan(plnParallel,scenarioId);
    scenarioDijs{scenarioIx} = matRad_calcDoseInfluence(ct,cst,stf,plnScenario);
end

dij = matRad_assembleParallelScenarioDijCore(scenarioDijs,scenarioIds, ...
    scenarioModel);
useParallel = true;

end

function scenarioModel = resolveScenarioModel(ct,pln,engine,matRad_cfg)
scenarioModel = [];

if isfield(pln,'multScen') && isa(pln.multScen,'matRad_ScenarioModel')
    scenarioModel = pln.multScen;
elseif isa(engine.multScen,'matRad_ScenarioModel')
    scenarioModel = engine.multScen;
else
    try
        scenarioModel = matRad_createScenarioModel(ct,engine.multScen);
    catch
        matRad_cfg.dispWarning(['UseParallel was requested for ', ...
            'multi-scenario dij calculation, but no valid scenario ', ...
            'model could be resolved. Falling back to serial dose ', ...
            'calculation.\n']);
    end
end
end

function workerMemoryBytes = estimateScenarioDoseWorkerMemoryBytes(ct,cst,stf,pln,engine)
inputBytes = matRad_variableBytes(ct) + matRad_variableBytes(cst) + ...
    matRad_variableBytes(stf) + matRad_variableBytes(pln);

numBixels = estimateNumBixels(stf);
numDoseVoxels = estimateNumDoseVoxels(ct,engine);
matrixBytes = double(numDoseVoxels) * double(max(1,numBixels)) * 8;

% Sparse dijs are usually far smaller than dense matrices, but each worker
% also holds temporary ray-tracing and sparse assembly state.
workerMemoryBytes = inputBytes + 0.05 * matrixBytes + ...
    max(128 * 1024^2,double(numDoseVoxels) * 8 * 10);
end

function numBixels = estimateNumBixels(stf)
if isfield(stf,'totalNumOfBixels')
    numBixels = sum([stf(:).totalNumOfBixels]);
else
    numBixels = 1;
end
end

function numDoseVoxels = estimateNumDoseVoxels(ct,engine)
doseGrid = engine.doseGrid;

if all(isfield(doseGrid,{'x','y','z'}))
    numDoseVoxels = numel(doseGrid.x) * numel(doseGrid.y) * ...
        numel(doseGrid.z);
    return;
end

ctGrid = matRad_getWorldAxes(ct);
numDoseVoxels = numel(ctGrid.x(1):doseGrid.resolution.x:ctGrid.x(end)) * ...
    numel(ctGrid.y(1):doseGrid.resolution.y:ctGrid.y(end)) * ...
    numel(ctGrid.z(1):doseGrid.resolution.z:ctGrid.z(end));
end
