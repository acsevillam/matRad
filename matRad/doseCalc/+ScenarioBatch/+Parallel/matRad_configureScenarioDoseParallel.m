function [useParallel, parallelProvider, parallelPlan] = matRad_configureScenarioDoseParallel( ...
                                                                                              provider, ctx, cfg, matRadCfg, stageName, ...
                                                                                              workerMemoryBytes, resultBytesPerScenario, ...
                                                                                              accumulatorBytes)
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
if nargin < 7 || isempty(resultBytesPerScenario)
    resultBytesPerScenario = 0;
end
if nargin < 8 || isempty(accumulatorBytes)
    accumulatorBytes = 0;
end

useParallel = false;
parallelProvider = provider;
parallelPlan = struct('useParallel', false, 'chunkSize', 1, ...
                      'workerUpperBound', 1, 'fallbackReason', '');

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
    workerMemoryBytes = ScenarioBatch.Parallel.matRad_estimateScenarioDoseWorkerMemoryBytes( ...
                                                                                            parallelProvider, provider, ctx);
end
parallelPlan = ScenarioBatch.Parallel.matRad_planScenarioDoseParallelStage( ...
                                                                           cfg, numScenarios, stageName, workerMemoryBytes, ...
                                                                           resultBytesPerScenario, accumulatorBytes, matRadCfg);
if ~parallelPlan.useParallel
    matRadCfg.dispWarning(['UseParallel was requested for %s, but the ', ...
                           'scenario-stage memory plan only allows one concurrent ', ...
                           'scenario (%s). Falling back to serial scenario-batch.\n'], ...
                          stageName, parallelPlan.fallbackReason);
    return
end

parallelOptions = ScenarioBatch.Pool.matRad_doseParallelPoolOptions( ...
                                                                    cfg, matRadCfg, 'parallelOptions');
[useParallel, ~, ~] = ScenarioBatch.Pool.matRad_configureSafeDoseParallelPool( ...
                                                                              workerMemoryBytes, numScenarios, matRadCfg, stageName, ...
                                                                              'fallbackDescription', 'serial scenario-batch', ...
                                                                              'memoryBudgetBytes', parallelPlan.memoryBudgetBytes, ...
                                                                              'workerUpperBound', parallelPlan.workerUpperBound, ...
                                                                              parallelOptions{:});
parallelPlan.useParallel = useParallel;
end

function provider = matRad_stripPreloadedScenarioDij(provider)
fieldsToRemove = {'preloadedScenarioId', 'preloadedDij'};
for i = 1:numel(fieldsToRemove)
    if isfield(provider, fieldsToRemove{i})
        provider = rmfield(provider, fieldsToRemove{i});
    end
end
end
