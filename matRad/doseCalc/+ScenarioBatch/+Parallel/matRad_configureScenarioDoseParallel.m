function [useParallel, parallelProvider] = matRad_configureScenarioDoseParallel(provider, ctx, cfg, matRadCfg, stageName, workerMemoryBytes)
% matRad_configureScenarioDoseParallel configures scenario-parallel precomputation
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

if nargin < 5 || isempty(stageName)
    stageName = 'scenario-batch scenario dose';
end
if nargin < 6
    workerMemoryBytes = [];
end

useParallel = false;
parallelProvider = provider;

if ~isfield(cfg, 'UseParallel') || ~cfg.UseParallel
    return
end

numScenarios = numel(ctx.scenarioDijIx);
if numScenarios < 2
    matRadCfg.dispWarning(['UseParallel was requested for %s, but fewer ', ...
                           'than two scenarios are selected. Falling back to serial ', ...
                           'scenario-batch.\n'], stageName);
    return
end

if isfield(provider, 'type') && strcmp(provider.type, 'inMemory')
    matRadCfg.dispWarning(['UseParallel was requested for %s, but ', ...
                           'precomputed dij input is used. Falling back to serial scenario-batch ', ...
                           'to avoid broadcasting large dose matrices to workers.\n'], stageName);
    return
end

parallelProvider = matRad_stripPreloadedScenarioDij(provider);
if isfield(parallelProvider, 'pln')
    parallelProvider.pln = ScenarioBatch.Worker.matRad_sanitizeWorkerPlan(parallelProvider.pln);
end
if isempty(workerMemoryBytes)
    workerMemoryBytes = matRad_estimateScenarioDoseWorkerMemoryBytes( ...
                                                                     parallelProvider, provider, ctx);
end

parallelOptions = ScenarioBatch.Pool.matRad_doseParallelPoolOptions( ...
                                                                    cfg, matRadCfg, 'parallelOptions');
[useParallel, ~, ~] = ScenarioBatch.Pool.matRad_configureSafeDoseParallelPool( ...
                                                                              workerMemoryBytes, numScenarios, matRadCfg, stageName, ...
                                                                              'fallbackDescription', 'serial scenario-batch', ...
                                                                              parallelOptions{:});
end

function provider = matRad_stripPreloadedScenarioDij(provider)
fieldsToRemove = {'preloadedScenarioId', 'preloadedDij'};
for i = 1:numel(fieldsToRemove)
    if isfield(provider, fieldsToRemove{i})
        provider = rmfield(provider, fieldsToRemove{i});
    end
end
end

function workerMemoryBytes = matRad_estimateScenarioDoseWorkerMemoryBytes( ...
                                                                          parallelProvider, originalProvider, ctx)
providerBytes = ScenarioBatch.Resources.matRad_estimateVariableBytes(parallelProvider);
scenarioDijBytes = 0;
if isfield(originalProvider, 'preloadedDij') && ~isempty(originalProvider.preloadedDij)
    scenarioDijBytes = ScenarioBatch.Resources.matRad_estimateVariableBytes(originalProvider.preloadedDij);
end

numRows = max([1; numel(ctx.targetRows); numel(ctx.oarRows)]);
% Scenario row workspaces keep sparse dose influence rows, matching the
% heuristic used for multi-scenario dij calculation.
sparseWorkspaceFillFactor = 0.05;
rowMatrixBytes = double(numRows) * double(max(1, ctx.numBixels)) * 8;
rowWorkspaceBytes = sparseWorkspaceFillFactor * rowMatrixBytes;
temporaryWorkspaceBytes = max(128 * 1024^2, double(numRows) * 8 * 10);
workerMemoryBytes = providerBytes + scenarioDijBytes + rowWorkspaceBytes + ...
    temporaryWorkspaceBytes;
end
