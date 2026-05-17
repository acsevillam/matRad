function [dij, useParallel] = matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine)
% matRad_calcDoseInfluenceParallelScenarios computes robust pencil-beam dijs in parallel
%
% call
%   [dij,useParallel] = ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct,cst,stf,pln,engine)
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

matRadCfg = MatRad_Config.instance();
dij = [];
useParallel = false;

if nargin < 5 || isempty(engine)
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
end

[parallelSupported, unsupportedReason] = ...
    ScenarioBatch.Parallel.matRad_supportsParallelScenarioDij(engine);
if ~parallelSupported
    switch unsupportedReason
        case 'notAnalyticalPencilBeam'
            matRadCfg.dispWarning(['UseParallel was requested for ', ...
                                   'multi-scenario dij calculation, but dose engine "%s" is ', ...
                                   'not an analytical pencil-beam engine. Falling back to ', ...
                                   'serial dose calculation.\n'], engine.name);
        case 'stochasticDijSampling'
            matRadCfg.dispWarning(['UseParallel was requested for ', ...
                                   'multi-scenario dij calculation, but dose engine "%s" ', ...
                                   'uses stochastic dij sampling. Falling back to serial ', ...
                                   'dose calculation.\n'], engine.name);
        otherwise
            matRadCfg.dispWarning(['UseParallel was requested for ', ...
                                   'multi-scenario dij calculation, but dose engine "%s" ', ...
                                   'does not support it. Falling back to serial dose ', ...
                                   'calculation.\n'], engine.name);
    end
    return
end

scenarioModel = matRad_resolveScenarioModel(ct, pln, engine, matRadCfg);
if isempty(scenarioModel)
    return
end

if scenarioModel.hasActiveAngularScenarioDimension() && ...
        scenarioModel.numOfBeams ~= numel(stf)
    scenarioModel.numOfBeams = numel(stf);
end
engine.multScen = scenarioModel;
engine.assertSupportedScenarioDimensions();

scenarioIds = scenarioModel.scenarioIds();
numScenarios = numel(scenarioIds);
if numScenarios < 2
    matRadCfg.dispWarning(['UseParallel was requested for multi-scenario ', ...
                           'dij calculation, but fewer than two scenarios are active. ', ...
                           'Falling back to serial dose calculation.\n']);
    return
end

plnParallel = ScenarioBatch.Worker.matRad_sanitizeWorkerPlan(pln, engine);
plnParallel.multScen = scenarioModel;

workerMemoryBytes = matRad_estimateScenarioDoseWorkerMemoryBytes(ct, cst, stf, ...
                                                                 plnParallel, engine);
parallelOptions = {};
if isfield(pln, 'propDoseCalc')
    parallelOptions = ScenarioBatch.Pool.matRad_doseParallelPoolOptions( ...
                                                                        pln.propDoseCalc, matRadCfg, 'propDoseCalc.parallelOptions');
end
[poolReady, ~, ~] = ScenarioBatch.Pool.matRad_configureSafeDoseParallelPool( ...
                                                                            workerMemoryBytes, numScenarios, matRadCfg, ...
                                                                            'multi-scenario dij calculation', ...
                                                                            'fallbackDescription', 'serial dose calculation', ...
                                                                            parallelOptions{:});
if ~poolReady
    return
end

scenarioDijs = cell(numScenarios, 1);

parfor scenarioIx = 1:numScenarios
    scenarioId = scenarioIds(scenarioIx);
    plnScenario = ScenarioBatch.Scenarios.matRad_buildSingleScenarioPlan(plnParallel, scenarioId);
    scenarioDijs{scenarioIx} = matRad_calcDoseInfluence(ct, cst, stf, plnScenario);
end

dij = ScenarioBatch.Dij.matRad_assembleParallelScenarioDij(scenarioDijs, scenarioIds, ...
                                                           scenarioModel);
useParallel = true;

end

function scenarioModel = matRad_resolveScenarioModel(ct, pln, engine, matRadCfg)
scenarioModel = [];

if isfield(pln, 'multScen') && isa(pln.multScen, 'matRad_ScenarioModel')
    scenarioModel = pln.multScen;
elseif isa(engine.multScen, 'matRad_ScenarioModel')
    scenarioModel = engine.multScen;
else
    try
        scenarioModel = matRad_ScenarioModel.create(engine.multScen, ct);
    catch
        matRadCfg.dispWarning(['UseParallel was requested for ', ...
                               'multi-scenario dij calculation, but no valid scenario ', ...
                               'model could be resolved. Falling back to serial dose ', ...
                               'calculation.\n']);
    end
end
end

function workerMemoryBytes = matRad_estimateScenarioDoseWorkerMemoryBytes(ct, cst, stf, pln, engine)
inputBytes = ScenarioBatch.Resources.matRad_estimateVariableBytes(ct) + ScenarioBatch.Resources.matRad_estimateVariableBytes(cst) + ...
    ScenarioBatch.Resources.matRad_estimateVariableBytes(stf) + ScenarioBatch.Resources.matRad_estimateVariableBytes(pln);

numBixels = matRad_estimateNumBixels(stf);
numDoseVoxels = matRad_estimateNumDoseVoxels(ct, engine);
matrixBytes = double(numDoseVoxels) * double(max(1, numBixels)) * 8;

% Sparse dijs are usually far smaller than dense matrices, but each worker
% also holds temporary ray-tracing and sparse assembly state.
workerMemoryBytes = inputBytes + 0.05 * matrixBytes + ...
    max(128 * 1024^2, double(numDoseVoxels) * 8 * 10);
end

function numBixels = matRad_estimateNumBixels(stf)
if isfield(stf, 'totalNumOfBixels')
    numBixels = sum([stf(:).totalNumOfBixels]);
else
    numBixels = 1;
end
end

function numDoseVoxels = matRad_estimateNumDoseVoxels(ct, engine)
doseGrid = engine.doseGrid;

if all(isfield(doseGrid, {'x', 'y', 'z'}))
    numDoseVoxels = numel(doseGrid.x) * numel(doseGrid.y) * ...
        numel(doseGrid.z);
    return
end

ctGrid = matRad_getWorldAxes(ct);
numDoseVoxels = numel(ctGrid.x(1):doseGrid.resolution.x:ctGrid.x(end)) * ...
    numel(ctGrid.y(1):doseGrid.resolution.y:ctGrid.y(end)) * ...
    numel(ctGrid.z(1):doseGrid.resolution.z:ctGrid.z(end));
end
