function parallelPlan = matRad_planScenarioDoseParallelStage( ...
                                                             cfg, numScenarios, stageName, workerBytes, ...
                                                             resultBytesPerScenario, accumulatorBytes, matRadCfg)
% matRad_planScenarioDoseParallelStage plans memory-safe scenario chunks
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

if nargin < 7 || isempty(matRadCfg)
    matRadCfg = MatRad_Config.instance();
end
if nargin < 3 || isempty(stageName)
    stageName = 'scenario-batch scenario dose';
end
if nargin < 4 || isempty(workerBytes)
    workerBytes = 0;
end
if nargin < 5 || isempty(resultBytesPerScenario)
    resultBytesPerScenario = 0;
end
if nargin < 6 || isempty(accumulatorBytes)
    accumulatorBytes = 0;
end

numScenarios = matRad_validatePositiveInteger(numScenarios, 'numScenarios', matRadCfg);
workerBytes = matRad_validateNonnegativeFiniteScalar(workerBytes, 'workerBytes', matRadCfg);
resultBytesPerScenario = matRad_validateNonnegativeFiniteScalar( ...
                                                                resultBytesPerScenario, ...
                                                                'resultBytesPerScenario', matRadCfg);
accumulatorBytes = matRad_validateNonnegativeFiniteScalar(accumulatorBytes, ...
                                                          'accumulatorBytes', matRadCfg);

[safetyFactor, minWorkerMemoryBytes, configuredWorkerUpperBound] = ...
    matRad_parsePlannerParallelOptions(cfg, matRadCfg);
effectiveWorkerBytes = max(workerBytes, minWorkerMemoryBytes) * safetyFactor;
memoryBudgetBytes = double(cfg.MemoryLimitMB) * 1e6;
taskBytes = effectiveWorkerBytes + resultBytesPerScenario;
availableForTasksBytes = memoryBudgetBytes - accumulatorBytes;

if taskBytes <= 0
    maxConcurrentByMemory = numScenarios;
elseif availableForTasksBytes <= 0
    maxConcurrentByMemory = 0;
else
    maxConcurrentByMemory = floor(availableForTasksBytes / taskBytes);
end
maxConcurrentByMemory = max(0, min(numScenarios, maxConcurrentByMemory));

[allocatedCpuCount, allocatedCpuSource] = matRad_getAllocatedCpuCount();
workerUpperBound = maxConcurrentByMemory;
if ~isempty(configuredWorkerUpperBound)
    workerUpperBound = min(workerUpperBound, configuredWorkerUpperBound);
end
if ~isempty(allocatedCpuCount)
    workerUpperBound = min(workerUpperBound, allocatedCpuCount);
end
workerUpperBound = min(workerUpperBound, numScenarios);

parallelPlan = struct();
parallelPlan.stageName = stageName;
parallelPlan.numScenarios = numScenarios;
parallelPlan.useParallel = workerUpperBound >= 2;
parallelPlan.chunkSize = max(1, min(numScenarios, workerUpperBound));
parallelPlan.workerUpperBound = max(1, workerUpperBound);
parallelPlan.fallbackReason = '';
parallelPlan.memoryBudgetBytes = memoryBudgetBytes;
parallelPlan.availableForTasksBytes = availableForTasksBytes;
parallelPlan.rawWorkerBytes = workerBytes;
parallelPlan.workerBytes = effectiveWorkerBytes;
parallelPlan.resultBytesPerScenario = resultBytesPerScenario;
parallelPlan.accumulatorBytes = accumulatorBytes;
parallelPlan.maxConcurrentByMemory = maxConcurrentByMemory;
parallelPlan.configuredWorkerUpperBound = configuredWorkerUpperBound;
parallelPlan.allocatedCpuCount = allocatedCpuCount;
parallelPlan.allocatedCpuSource = allocatedCpuSource;

if ~parallelPlan.useParallel
    if maxConcurrentByMemory < 2
        parallelPlan.fallbackReason = 'memoryBudget';
    elseif ~isempty(configuredWorkerUpperBound) && configuredWorkerUpperBound < 2
        parallelPlan.fallbackReason = 'workerUpperBound';
    elseif ~isempty(allocatedCpuCount) && allocatedCpuCount < 2
        parallelPlan.fallbackReason = allocatedCpuSource;
    else
        parallelPlan.fallbackReason = 'workerLimit';
    end
    return
end

matRadCfg.dispInfo(['matRad: %s scenario parallel plan uses chunks of ', ...
                    '%d scenario(s), worker cap %d, and memory budget %.2f MB.\n'], ...
                   stageName, parallelPlan.chunkSize, ...
                   parallelPlan.workerUpperBound, memoryBudgetBytes / 1e6);
matRadCfg.dispInfo(['matRad: %s scenario parallel memory model: worker %.2f MB, ', ...
                    'result %.2f MB/scenario, accumulator %.2f MB.\n'], ...
                   stageName, effectiveWorkerBytes / 1e6, ...
                   resultBytesPerScenario / 1e6, accumulatorBytes / 1e6);
end

function value = matRad_validatePositiveInteger(value, valueName, matRadCfg)
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
     round(value) == value && value >= 1)
    matRadCfg.dispError('%s must be a positive integer scalar.', valueName);
end
value = double(value);
end

function value = matRad_validateNonnegativeFiniteScalar(value, valueName, matRadCfg)
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0)
    matRadCfg.dispError('%s must be a nonnegative finite scalar.', valueName);
end
value = double(value);
end

function [safetyFactor, minWorkerMemoryBytes, workerUpperBound] = ...
    matRad_parsePlannerParallelOptions(cfg, matRadCfg)
safetyFactor = 1.2;
minWorkerMemoryBytes = 512 * 1024^2;
workerUpperBound = [];

if ~isfield(cfg, 'parallelOptions') || isempty(cfg.parallelOptions)
    return
end

parallelOptions = cfg.parallelOptions;
if ~isstruct(parallelOptions) || ~isscalar(parallelOptions)
    matRadCfg.dispError('parallelOptions must be a scalar struct.');
end

fields = fieldnames(parallelOptions);
for i = 1:numel(fields)
    fieldName = fields{i};
    value = parallelOptions.(fieldName);
    switch fieldName
        case 'workerMemorySafetyFactor'
            safetyFactor = matRad_validatePositiveFiniteScalar(value, ...
                                                               fieldName, matRadCfg);
        case 'minWorkerMemoryBytes'
            minWorkerMemoryBytes = matRad_validateNonnegativeFiniteScalar(value, ...
                                                                          fieldName, matRadCfg);
        case 'workerUpperBound'
            workerUpperBound = matRad_validateOptionalPositiveInteger(value, ...
                                                                      fieldName, matRadCfg);
        case 'memoryReserveFraction'
        otherwise
            matRadCfg.dispError('Unknown dose parallel option "%s" in parallelOptions.', ...
                                fieldName);
    end
end
end

function value = matRad_validatePositiveFiniteScalar(value, valueName, matRadCfg)
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 1)
    matRadCfg.dispError('%s must be a positive finite scalar.', valueName);
end
value = double(value);
end

function value = matRad_validateOptionalPositiveInteger(value, valueName, matRadCfg)
if isempty(value)
    return
end
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
     round(value) == value && value >= 1)
    matRadCfg.dispError('%s must be a positive integer scalar or empty.', ...
                        valueName);
end
value = double(value);
end
