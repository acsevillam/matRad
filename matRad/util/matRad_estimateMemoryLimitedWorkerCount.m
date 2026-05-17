function [maxWorkers, memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount(workerMemoryBytes, varargin)
% matRad_estimateMemoryLimitedWorkerCount estimates a memory-safe worker count
%
% call
%   [maxWorkers,memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount(workerMemoryBytes)
%   [maxWorkers,memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount(workerMemoryBytes,'numTasks',numTasks)
%   [maxWorkers,memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount(...,'safetyFactor',safetyFactor)
%   [maxWorkers,memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount(...,'reserveFraction',reserveFraction)
%   [maxWorkers,memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount(...,'workerUpperBound',workerUpperBound)
%   [maxWorkers,memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount(...,'minWorkerMemoryBytes',minWorkerMemoryBytes)
%
% input
%   workerMemoryBytes:     estimated memory footprint per worker in bytes
%
% input (optional Name-Value pairs)
%   varargin:             optional Name-Value pairs
%   numTasks:              maximum useful worker count imposed by the number
%                          of independent tasks
%   safetyFactor:          factor applied to workerMemoryBytes
%   reserveFraction:       fraction of total system memory kept in reserve
%   workerUpperBound:      explicit upper bound for the returned worker count
%   minWorkerMemoryBytes:  lower bound applied before safetyFactor
%   limitToDefaultPool:    keep the result at or below the default local pool
%                          size if workerUpperBound is empty
%
% output
%   maxWorkers:            maximum recommended worker count
%   memoryEstimate:        struct with diagnostic memory and worker estimates
%                          in bytes
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

p = matRad_parseMemoryLimitedWorkerInputs(workerMemoryBytes, varargin{:});

workerUpperBound = matRad_validateOptionalPositiveInteger( ...
                                                          p.Results.workerUpperBound, 'workerUpperBound');

maxWorkers = [];
memoryInfo = matRad_getSystemMemoryInfo('reserveFraction', p.Results.reserveFraction, ...
                                        'minReserveBytes', 1 * 1024^3);
memoryEstimate = memoryInfo;
memoryEstimate.rawWorkerBytes = workerMemoryBytes;
memoryEstimate.workerBytes = [];
memoryEstimate.memoryLimitedWorkers = [];
memoryEstimate.defaultPoolSize = [];
memoryEstimate.workerUpperBound = workerUpperBound;
memoryEstimate.numTasks = p.Results.numTasks;
memoryEstimate.safetyFactor = p.Results.safetyFactor;
memoryEstimate.reserveFraction = p.Results.reserveFraction;
memoryEstimate.minWorkerMemoryBytes = p.Results.minWorkerMemoryBytes;
memoryEstimate.limitToDefaultPool = logical(p.Results.limitToDefaultPool);
memoryEstimate.maxWorkers = [];

memoryEstimate.workerBytes = max(workerMemoryBytes, p.Results.minWorkerMemoryBytes) * p.Results.safetyFactor;
memoryEstimate.defaultPoolSize = matRad_getDefaultParallelPoolSize();

if isempty(memoryEstimate.usableBytes) || isempty(memoryEstimate.workerBytes) || memoryEstimate.workerBytes <= 0
    return
end

memoryEstimate.memoryLimitedWorkers = max(1, floor(memoryEstimate.usableBytes / memoryEstimate.workerBytes));

maxWorkers = memoryEstimate.memoryLimitedWorkers;

if ~isempty(workerUpperBound)
    maxWorkers = min(maxWorkers, workerUpperBound);
elseif p.Results.limitToDefaultPool && ~isempty(memoryEstimate.defaultPoolSize)
    maxWorkers = min(maxWorkers, memoryEstimate.defaultPoolSize);
end

if ~isempty(p.Results.numTasks) && isfinite(p.Results.numTasks)
    maxWorkers = min(maxWorkers, p.Results.numTasks);
end

memoryEstimate.maxWorkers = maxWorkers;
end

function p = matRad_parseMemoryLimitedWorkerInputs(workerMemoryBytes, varargin)
p = inputParser;
p.addRequired('workerMemoryBytes', @matRad_isNonnegativeFiniteScalar);
p.addParameter('numTasks', Inf, @matRad_isEmptyOrPositiveScalar);
p.addParameter('safetyFactor', 1, @matRad_isPositiveFiniteScalar);
p.addParameter('reserveFraction', 0.10, @matRad_isReserveFraction);
p.addParameter('workerUpperBound', [], @(x) true);
p.addParameter('minWorkerMemoryBytes', 0, @matRad_isNonnegativeFiniteScalar);
p.addParameter('limitToDefaultPool', true, @matRad_isScalarLogicalLike);
p.parse(workerMemoryBytes, varargin{:});
end

function value = matRad_validateOptionalPositiveInteger(value, valueName)
if isempty(value)
    return
end

if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
     round(value) == value && value >= 1)
    matRadCfg = MatRad_Config.instance();
    matRadCfg.dispError('%s must be a positive integer scalar or empty.', ...
                        valueName);
end
end

function tf = matRad_isNonnegativeFiniteScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0;
end

function tf = matRad_isPositiveFiniteScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && value >= 1;
end

function tf = matRad_isEmptyOrPositiveScalar(value)
tf = isempty(value) || (isnumeric(value) && isscalar(value) && value >= 1);
end

function tf = matRad_isReserveFraction(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0 && value < 1;
end

function tf = matRad_isScalarLogicalLike(value)
tf = (islogical(value) || isnumeric(value)) && isscalar(value);
end

function defaultPoolSize = matRad_getDefaultParallelPoolSize()
defaultPoolSize = [];
try
    cluster = parcluster();
    defaultPoolSize = cluster.NumWorkers;
catch
end

if isempty(defaultPoolSize) || ~isnumeric(defaultPoolSize) || defaultPoolSize < 1
    try
        defaultPoolSize = feature('numcores');
    catch
        defaultPoolSize = [];
    end
end
end
