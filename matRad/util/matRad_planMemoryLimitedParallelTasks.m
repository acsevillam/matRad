function parallelPlan = matRad_planMemoryLimitedParallelTasks(numTasks, workerBytes, varargin)
% matRad_planMemoryLimitedParallelTasks plans memory-safe parallel task chunks
%
% call
%   parallelPlan = matRad_planMemoryLimitedParallelTasks(numTasks,workerBytes)
%   parallelPlan = matRad_planMemoryLimitedParallelTasks(...,'resultBytesPerTask',resultBytesPerTask)
%   parallelPlan = matRad_planMemoryLimitedParallelTasks(...,'accumulatorBytes',accumulatorBytes)
%   parallelPlan = matRad_planMemoryLimitedParallelTasks(...,'memoryBudgetBytes',memoryBudgetBytes)
%
% input
%   numTasks:    number of independent tasks
%   workerBytes: estimated memory footprint per worker in bytes
%
% output
%   parallelPlan: struct with worker and chunk limits
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

p = matRad_parseMemoryLimitedParallelTaskInputs(varargin{:});
matRadCfg = p.Results.matRadCfg;
if isempty(matRadCfg)
    matRadCfg = MatRad_Config.instance();
end

numTasks = matRad_validatePositiveInteger(numTasks, 'numTasks', matRadCfg);
workerBytes = matRad_validateNonnegativeFiniteScalar(workerBytes, ...
                                                    'workerBytes', matRadCfg);
resultBytesPerTask = matRad_validateNonnegativeFiniteScalar( ...
    p.Results.resultBytesPerTask, 'resultBytesPerTask', matRadCfg);
accumulatorBytes = matRad_validateNonnegativeFiniteScalar( ...
    p.Results.accumulatorBytes, 'accumulatorBytes', matRadCfg);
safetyFactor = matRad_validatePositiveFiniteScalar( ...
    p.Results.safetyFactor, 'safetyFactor', matRadCfg);
minWorkerMemoryBytes = matRad_validateNonnegativeFiniteScalar( ...
    p.Results.minWorkerMemoryBytes, 'minWorkerMemoryBytes', matRadCfg);
configuredWorkerUpperBound = matRad_validateOptionalPositiveInteger( ...
    p.Results.workerUpperBound, 'workerUpperBound', matRadCfg);

effectiveWorkerBytes = max(workerBytes, minWorkerMemoryBytes) * safetyFactor;
[memoryBudgetBytes, memoryInfo] = matRad_resolveParallelTaskMemoryBudget( ...
    p.Results.memoryBudgetBytes, p.Results.memoryInfo, p.Results.reserveFraction, ...
    p.Results.minReserveBytes, p.Results.environment, p.Results.cgroupRoot, matRadCfg);

taskBytes = effectiveWorkerBytes + resultBytesPerTask;
if isempty(memoryBudgetBytes) || ~isfinite(memoryBudgetBytes) || memoryBudgetBytes <= 0
    availableForTasksBytes = [];
    maxConcurrentByMemory = 0;
else
    availableForTasksBytes = memoryBudgetBytes - accumulatorBytes;
    if taskBytes <= 0
        maxConcurrentByMemory = numTasks;
    elseif availableForTasksBytes <= 0
        maxConcurrentByMemory = 0;
    else
        maxConcurrentByMemory = floor(availableForTasksBytes / taskBytes);
    end
end
maxConcurrentByMemory = max(0, min(numTasks, maxConcurrentByMemory));

[allocatedCpuCount, allocatedCpuSource] = matRad_getAllocatedCpuCount( ...
    'environment', p.Results.environment);

workerUpperBound = maxConcurrentByMemory;
if ~isempty(configuredWorkerUpperBound)
    workerUpperBound = min(workerUpperBound, configuredWorkerUpperBound);
end
if ~isempty(allocatedCpuCount)
    workerUpperBound = min(workerUpperBound, allocatedCpuCount);
end
workerUpperBound = min(workerUpperBound, numTasks);

parallelPlan = struct();
parallelPlan.stageName = char(p.Results.stageName);
parallelPlan.numTasks = numTasks;
parallelPlan.useParallel = workerUpperBound >= 2;
parallelPlan.chunkSize = max(1, min(numTasks, workerUpperBound));
parallelPlan.workerUpperBound = max(1, workerUpperBound);
parallelPlan.fallbackReason = '';
parallelPlan.memoryInfo = memoryInfo;
parallelPlan.memoryBudgetBytes = memoryBudgetBytes;
parallelPlan.availableForTasksBytes = availableForTasksBytes;
parallelPlan.rawWorkerBytes = workerBytes;
parallelPlan.effectiveWorkerBytes = effectiveWorkerBytes;
parallelPlan.workerBytes = effectiveWorkerBytes;
parallelPlan.resultBytesPerTask = resultBytesPerTask;
parallelPlan.resultBytesPerScenario = resultBytesPerTask;
parallelPlan.accumulatorBytes = accumulatorBytes;
parallelPlan.maxConcurrentByMemory = maxConcurrentByMemory;
parallelPlan.configuredWorkerUpperBound = configuredWorkerUpperBound;
parallelPlan.allocatedCpuCount = allocatedCpuCount;
parallelPlan.allocatedCpuSource = allocatedCpuSource;
parallelPlan.safetyFactor = safetyFactor;
parallelPlan.minWorkerMemoryBytes = minWorkerMemoryBytes;
parallelPlan.reserveFraction = p.Results.reserveFraction;

if ~parallelPlan.useParallel
    parallelPlan.fallbackReason = matRad_parallelTaskFallbackReason( ...
        memoryBudgetBytes, maxConcurrentByMemory, configuredWorkerUpperBound, ...
        allocatedCpuCount, allocatedCpuSource);
end
end

function p = matRad_parseMemoryLimitedParallelTaskInputs(varargin)
p = inputParser;
p.addParameter('stageName', 'parallel tasks', ...
               @(x) ischar(x) || (isstring(x) && isscalar(x)));
p.addParameter('resultBytesPerTask', 0, @matRad_isNonnegativeFiniteScalar);
p.addParameter('accumulatorBytes', 0, @matRad_isNonnegativeFiniteScalar);
p.addParameter('memoryBudgetBytes', [], @matRad_isEmptyOrPositiveFiniteScalar);
p.addParameter('memoryInfo', [], @(x) isempty(x) || isstruct(x));
p.addParameter('reserveFraction', 0.10, @matRad_isReserveFraction);
p.addParameter('minReserveBytes', 1 * 1024^3, @matRad_isNonnegativeFiniteScalar);
p.addParameter('safetyFactor', 1, @matRad_isPositiveFiniteScalar);
p.addParameter('minWorkerMemoryBytes', 0, @matRad_isNonnegativeFiniteScalar);
p.addParameter('workerUpperBound', [], @(x) true);
p.addParameter('environment', [], @(x) isempty(x) || isstruct(x));
p.addParameter('cgroupRoot', '/sys/fs/cgroup', ...
               @(x) ischar(x) || (isstring(x) && isscalar(x)));
p.addParameter('matRadCfg', [], @(x) isempty(x) || isa(x, 'MatRad_Config'));
p.parse(varargin{:});
end

function [memoryBudgetBytes, memoryInfo] = matRad_resolveParallelTaskMemoryBudget( ...
    memoryBudgetBytes, memoryInfo, reserveFraction, minReserveBytes, environment, ...
    cgroupRoot, matRadCfg)
if ~isempty(memoryBudgetBytes)
    memoryBudgetBytes = matRad_validatePositiveFiniteScalar(memoryBudgetBytes, ...
                                                            'memoryBudgetBytes', matRadCfg);
    memoryInfo = struct('availableBytes', memoryBudgetBytes, ...
                        'totalBytes', memoryBudgetBytes, ...
                        'reserveBytes', 0, ...
                        'usableBytes', memoryBudgetBytes, ...
                        'source', 'explicit:memoryBudgetBytes');
    return
end

if isempty(memoryInfo)
    memoryInfo = matRad_getSystemMemoryInfo('reserveFraction', reserveFraction, ...
                                            'minReserveBytes', minReserveBytes, ...
                                            'environment', environment, ...
                                            'cgroupRoot', cgroupRoot);
end

memoryBudgetBytes = [];
if isstruct(memoryInfo) && isfield(memoryInfo, 'usableBytes') && ...
        ~isempty(memoryInfo.usableBytes) && isfinite(memoryInfo.usableBytes) && ...
        memoryInfo.usableBytes > 0
    memoryBudgetBytes = double(memoryInfo.usableBytes);
end
end

function reason = matRad_parallelTaskFallbackReason(memoryBudgetBytes, ...
                                                    maxConcurrentByMemory, ...
                                                    configuredWorkerUpperBound, ...
                                                    allocatedCpuCount, allocatedCpuSource)
if isempty(memoryBudgetBytes) || ~isfinite(memoryBudgetBytes) || memoryBudgetBytes <= 0
    reason = 'memoryBudgetUnavailable';
elseif maxConcurrentByMemory < 2
    reason = 'memoryBudget';
elseif ~isempty(configuredWorkerUpperBound) && configuredWorkerUpperBound < 2
    reason = 'workerUpperBound';
elseif ~isempty(allocatedCpuCount) && allocatedCpuCount < 2
    reason = allocatedCpuSource;
else
    reason = 'workerLimit';
end
end

function value = matRad_validatePositiveInteger(value, valueName, matRadCfg)
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
     round(value) == value && value >= 1)
    matRadCfg.dispError('%s must be a positive integer scalar.', valueName);
end
value = double(value);
end

function value = matRad_validateNonnegativeFiniteScalar(value, valueName, matRadCfg)
if ~matRad_isNonnegativeFiniteScalar(value)
    matRadCfg.dispError('%s must be a nonnegative finite scalar.', valueName);
end
value = double(value);
end

function value = matRad_validatePositiveFiniteScalar(value, valueName, matRadCfg)
if ~matRad_isPositiveFiniteScalar(value)
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

function tf = matRad_isNonnegativeFiniteScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0;
end

function tf = matRad_isPositiveFiniteScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && value > 0;
end

function tf = matRad_isEmptyOrPositiveFiniteScalar(value)
tf = isempty(value) || matRad_isPositiveFiniteScalar(value);
end

function tf = matRad_isReserveFraction(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0 && value < 1;
end
