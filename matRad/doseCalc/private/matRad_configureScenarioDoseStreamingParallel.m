function [useParallel,parallelProvider] = matRad_configureScenarioDoseStreamingParallel(provider,ctx,cfg,matRad_cfg,stageName,workerMemoryBytes)
% matRad_configureScenarioDoseStreamingParallel configures scenario-parallel streaming
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
    stageName = 'streaming scenario dose';
end
if nargin < 6
    workerMemoryBytes = [];
end

useParallel = false;
parallelProvider = provider;

if ~isfield(cfg,'UseParallel') || ~cfg.UseParallel
    return;
end

numScenarios = numel(ctx.scenarioDijIx);
if numScenarios < 2
    matRad_cfg.dispWarning(['UseParallel was requested for %s, but fewer ', ...
        'than two scenarios are selected. Falling back to serial ', ...
        'streaming.\n'],stageName);
    return;
end

if isfield(provider,'type') && strcmp(provider.type,'inMemory')
    matRad_cfg.dispWarning(['UseParallel was requested for %s, but ', ...
        'precomputed dij input is used. Falling back to serial streaming ', ...
        'to avoid broadcasting large dose matrices to workers.\n'],stageName);
    return;
end

parallelProvider = stripPreloadedScenarioDij(provider);
if isfield(parallelProvider,'pln')
    parallelProvider.pln = matRad_makeWorkerSafePlan(parallelProvider.pln);
end
if isempty(workerMemoryBytes)
    workerMemoryBytes = estimateStreamingScenarioWorkerMemoryBytes( ...
        parallelProvider,provider,ctx);
end

[useParallel,~,~] = matRad_configureSafeDoseParallelPool( ...
    workerMemoryBytes,numScenarios,matRad_cfg,stageName, ...
    'fallbackDescription','serial streaming');
end

function provider = stripPreloadedScenarioDij(provider)
fieldsToRemove = {'preloadedScenarioId','preloadedDij'};
for i = 1:numel(fieldsToRemove)
    if isfield(provider,fieldsToRemove{i})
        provider = rmfield(provider,fieldsToRemove{i});
    end
end
end

function workerMemoryBytes = estimateStreamingScenarioWorkerMemoryBytes( ...
    parallelProvider,originalProvider,ctx)
providerBytes = matRad_variableBytes(parallelProvider);
scenarioDijBytes = 0;
if isfield(originalProvider,'preloadedDij') && ~isempty(originalProvider.preloadedDij)
    scenarioDijBytes = matRad_variableBytes(originalProvider.preloadedDij);
end

numRows = max([1; numel(ctx.targetRows); numel(ctx.oarRows)]);
% Streaming row workspaces keep sparse dose influence rows, matching the
% heuristic used for multi-scenario dij calculation.
sparseWorkspaceFillFactor = 0.05;
rowMatrixBytes = double(numRows) * double(max(1,ctx.numBixels)) * 8;
rowWorkspaceBytes = sparseWorkspaceFillFactor * rowMatrixBytes;
temporaryWorkspaceBytes = max(128 * 1024^2,double(numRows) * 8 * 10);
workerMemoryBytes = providerBytes + scenarioDijBytes + rowWorkspaceBytes + ...
    temporaryWorkspaceBytes;
end
