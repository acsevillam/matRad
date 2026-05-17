function [useParallel, pPool, memoryEstimate] = matRad_configureSafeDoseParallelPool( ...
                                                                                     workerMemoryBytes, numTasks, matRadCfg, ...
                                                                                     calculationName, varargin)
% matRad_configureSafeDoseParallelPool configures safe dose parallelism
%   The helper may create, reduce, or increase the active parallel pool to
%   match the memory-safe worker count for the requested dose calculation.
%
% call
%   [useParallel,pPool,memoryEstimate] = ScenarioBatch.Pool.matRad_configureSafeDoseParallelPool( ...
%       workerMemoryBytes,numTasks,matRadCfg,calculationName)
%
% input
%   workerMemoryBytes: estimated bytes needed by one worker
%   numTasks:          number of independent parallel tasks
%   matRadCfg:        MatRad_Config instance
%   calculationName:   text used in diagnostics
%   varargin:          optional Name-Value pair arguments
%
% input (optional Name-Value pairs)
%   fallbackDescription:  text used when serial execution is selected
%   safetyFactor:         factor applied to the per-worker estimate
%   reserveFraction:      fraction of system memory kept in reserve
%   minWorkerMemoryBytes: lower bound applied before safetyFactor
%   workerUpperBound:     explicit upper bound for the worker count
%
% output
%   useParallel:       true if a pool with at least two workers is ready
%   pPool:             configured parallel pool, or []
%   memoryEstimate:    memory estimate returned by
%                      matRad_estimateMemoryLimitedWorkerCount
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

if nargin < 4 || isempty(calculationName)
    calculationName = 'dose calculation';
end

options = matRad_parseOptions(varargin{:});
useParallel = false;
pPool = [];
memoryEstimate = struct('workerBytes', []);

if numTasks < 2
    matRadCfg.dispWarning(['UseParallel was requested for %s, but fewer ', ...
                           'than two tasks are available. Falling back to %s.\n'], ...
                          calculationName, options.fallbackDescription);
    return
end

if ~matRad_isParallelComputingAvailable()
    matRadCfg.dispWarning(['UseParallel was requested for %s, but the ', ...
                           'Parallel Computing Toolbox is unavailable. Falling back to %s.\n'], ...
                          calculationName, options.fallbackDescription);
    return
end

[poolSizeLimit, memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount( ...
                                                                          workerMemoryBytes, 'numTasks', numTasks, ...
                                                                          'safetyFactor', options.safetyFactor, ...
                                                                          'reserveFraction', options.reserveFraction, ...
                                                                          'minWorkerMemoryBytes', options.minWorkerMemoryBytes, ...
                                                                          'workerUpperBound', options.workerUpperBound);

if isfield(memoryEstimate, 'workerBytes') && ...
        ~isempty(memoryEstimate.workerBytes)
    matRadCfg.dispInfo(['matRad: Estimated %s memory per parallel ', ...
                        'worker is %.2f MB.\n'], calculationName, ...
                       memoryEstimate.workerBytes / 1e6);
end

if isempty(poolSizeLimit) || poolSizeLimit < 2
    matRadCfg.dispWarning(['UseParallel was requested for %s, but the ', ...
                           'estimated memory only allows one worker. Falling back to %s.\n'], ...
                          calculationName, options.fallbackDescription);
    return
end

try
    pPool = matRad_configureParallelPoolSize(poolSizeLimit, ...
                                             calculationName, matRadCfg);
catch
    matRadCfg.dispWarning(['Could not configure a parallel pool for %s. ', ...
                           'Falling back to %s.\n'], calculationName, ...
                          options.fallbackDescription);
    return
end

if isempty(pPool) || pPool.NumWorkers < 2
    matRadCfg.dispWarning(['UseParallel was requested for %s, but fewer ', ...
                           'than two workers are available. Falling back to %s.\n'], ...
                          calculationName, options.fallbackDescription);
    return
end

matRadCfg.dispInfo(['matRad: Parallel %s uses %d worker(s) for %d ', ...
                    'task(s).\n'], calculationName, pPool.NumWorkers, numTasks);
useParallel = true;
end

function options = matRad_parseOptions(varargin)
options = struct();
options.fallbackDescription = 'serial dose calculation';
options.safetyFactor = 1.2;
options.reserveFraction = 0.10;
options.minWorkerMemoryBytes = 512 * 1024^2;
options.workerUpperBound = [];

if mod(numel(varargin), 2) ~= 0
    matRadCfg = MatRad_Config.instance();
    matRadCfg.dispError('Parallel pool options must be name-value pairs.');
end

for i = 1:2:numel(varargin)
    name = varargin{i};
    value = varargin{i + 1};
    if isstring(name) && isscalar(name)
        name = char(name);
    end
    switch lower(name)
        case 'fallbackdescription'
            options.fallbackDescription = value;
        case 'safetyfactor'
            options.safetyFactor = value;
        case 'reservefraction'
            options.reserveFraction = value;
        case 'minworkermemorybytes'
            options.minWorkerMemoryBytes = value;
        case 'workerupperbound'
            matRad_validateOptionalPositiveInteger(value, 'workerUpperBound');
            options.workerUpperBound = value;
        otherwise
            matRadCfg = MatRad_Config.instance();
            matRadCfg.dispError('Unknown parallel pool option "%s".', name);
    end
end
end

function matRad_validateOptionalPositiveInteger(value, valueName)
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

function available = matRad_isParallelComputingAvailable()
available = false;
if exist('parpool', 'file') ~= 2 || exist('gcp', 'file') ~= 2
    return
end

try
    [available, ~] = license('checkout', 'Distrib_Computing_Toolbox');
catch
    available = false;
end

if isempty(available)
    available = false;
end
available = logical(available);
end
