function pPool = matRad_configureParallelPoolSize(targetWorkers, ...
                                                  calculationName, matRadCfg)
% matRad_configureParallelPoolSize configures a parallel pool target size
%
% call
%   pPool = matRad_configureParallelPoolSize(targetWorkers,calculationName, ...
%       matRadCfg)
%
% input
%   targetWorkers:   positive integer number of workers requested by the
%                    caller after its own safety checks
%   calculationName: text used in diagnostics
%   matRadCfg:      MatRad_Config instance
%
% output
%   pPool:           configured parallel pool
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

if nargin < 2 || isempty(calculationName)
    calculationName = 'parallel calculation';
end
if nargin < 3 || isempty(matRadCfg)
    matRadCfg = MatRad_Config.instance();
end

matRad_validateTargetWorkers(targetWorkers, matRadCfg);
targetWorkers = double(targetWorkers);

pPool = gcp('nocreate');

if isempty(pPool)
    try
        pPool = parpool(targetWorkers);
    catch ME
        matRadCfg.dispError(['Could not create a %d-worker parallel ', ...
                             'pool for %s: %s'], targetWorkers, calculationName, ME.message);
    end
    return
end

currentWorkers = pPool.NumWorkers;
if currentWorkers == targetWorkers
    return
end

if currentWorkers > targetWorkers
    matRadCfg.dispWarning(['matRad is reducing the active parallel pool ', ...
                           'from %d to %d worker(s) for %s.\n'], currentWorkers, ...
                          targetWorkers, calculationName);
else
    matRadCfg.dispInfo(['matRad: Increasing the active parallel pool ', ...
                        'from %d to %d worker(s) for %s.\n'], currentWorkers, ...
                       targetWorkers, calculationName);
end

try
    delete(pPool);
    pPool = parpool(targetWorkers);
catch ME
    matRad_restoreOriginalPool(currentWorkers, targetWorkers, calculationName, ...
                               matRadCfg, ME);
end
end

function matRad_validateTargetWorkers(targetWorkers, matRadCfg)
if ~(isnumeric(targetWorkers) && isscalar(targetWorkers) && ...
     isfinite(targetWorkers) && round(targetWorkers) == targetWorkers && ...
     targetWorkers >= 1)
    matRadCfg.dispError('targetWorkers must be a positive integer scalar.');
end
end

function matRad_restoreOriginalPool(originalWorkers, targetWorkers, calculationName, ...
                                    matRadCfg, resizeError)
try
    pPool = gcp('nocreate');
    if ~isempty(pPool)
        delete(pPool);
    end
    parpool(originalWorkers);
catch restoreError
    matRadCfg.dispError(['Could not configure a %d-worker parallel pool ', ...
                         'for %s and could not restore the previous %d-worker pool. ', ...
                         'Resize error: %s Restore error: %s'], targetWorkers, calculationName, ...
                        originalWorkers, resizeError.message, restoreError.message);
end

matRadCfg.dispError(['Could not configure a %d-worker parallel pool for ', ...
                     '%s. The previous %d-worker pool was restored. Resize error: %s'], ...
                    targetWorkers, calculationName, originalWorkers, resizeError.message);
end
