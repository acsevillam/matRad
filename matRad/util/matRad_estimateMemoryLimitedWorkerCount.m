function [maxWorkers,memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount(workerMemoryBytes,varargin)
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

p = inputParser;
p.addRequired('workerMemoryBytes',@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('numTasks',Inf,@(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
p.addParameter('safetyFactor',1,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
p.addParameter('reserveFraction',0.10,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x < 1);
p.addParameter('workerUpperBound',[],@(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
p.addParameter('minWorkerMemoryBytes',0,@(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('limitToDefaultPool',true,@(x) (islogical(x) || isnumeric(x)) && isscalar(x));
p.parse(workerMemoryBytes,varargin{:});

maxWorkers = [];
memoryInfo = matRad_getSystemMemoryInfo('reserveFraction',p.Results.reserveFraction, ...
    'minReserveBytes',1 * 1024^3);
memoryEstimate = memoryInfo;
memoryEstimate.rawWorkerBytes = workerMemoryBytes;
memoryEstimate.workerBytes = [];
memoryEstimate.memoryLimitedWorkers = [];
memoryEstimate.defaultPoolSize = [];
memoryEstimate.workerUpperBound = p.Results.workerUpperBound;
memoryEstimate.numTasks = p.Results.numTasks;
memoryEstimate.safetyFactor = p.Results.safetyFactor;
memoryEstimate.reserveFraction = p.Results.reserveFraction;
memoryEstimate.minWorkerMemoryBytes = p.Results.minWorkerMemoryBytes;
memoryEstimate.limitToDefaultPool = logical(p.Results.limitToDefaultPool);
memoryEstimate.maxWorkers = [];

memoryEstimate.workerBytes = max(workerMemoryBytes,p.Results.minWorkerMemoryBytes) * p.Results.safetyFactor;
memoryEstimate.defaultPoolSize = getDefaultParallelPoolSize();

if isempty(memoryEstimate.usableBytes) || isempty(memoryEstimate.workerBytes) || memoryEstimate.workerBytes <= 0
    return;
end

memoryEstimate.memoryLimitedWorkers = max(1,floor(memoryEstimate.usableBytes / memoryEstimate.workerBytes));

maxWorkers = memoryEstimate.memoryLimitedWorkers;

if ~isempty(p.Results.workerUpperBound)
    maxWorkers = min(maxWorkers,p.Results.workerUpperBound);
elseif p.Results.limitToDefaultPool && ~isempty(memoryEstimate.defaultPoolSize)
    maxWorkers = min(maxWorkers,memoryEstimate.defaultPoolSize);
end

if ~isempty(p.Results.numTasks) && isfinite(p.Results.numTasks)
    maxWorkers = min(maxWorkers,p.Results.numTasks);
end

memoryEstimate.maxWorkers = maxWorkers;
end

function defaultPoolSize = getDefaultParallelPoolSize()
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
