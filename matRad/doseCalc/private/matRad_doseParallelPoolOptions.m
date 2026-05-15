function options = matRad_doseParallelPoolOptions(owner,matRad_cfg,context)
% matRad_doseParallelPoolOptions maps resource contracts to pool options
%
% call
%   options = matRad_doseParallelPoolOptions(owner,matRad_cfg,context)
%
% input
%   owner:      struct or object with a parallelOptions member
%   matRad_cfg: MatRad_Config instance for diagnostics
%   context:    text used in diagnostics
%
% output
%   options:    name-value cell array for matRad_configureSafeDoseParallelPool
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

if nargin < 2 || isempty(matRad_cfg)
    matRad_cfg = MatRad_Config.instance();
end
if nargin < 3 || isempty(context)
    context = 'parallelOptions';
end

options = {};
parallelOptions = [];
if isempty(owner)
    return;
elseif isstruct(owner)
    if ~isfield(owner,'parallelOptions')
        return;
    end
    parallelOptions = owner.parallelOptions;
elseif isobject(owner) && isprop(owner,'parallelOptions')
    parallelOptions = owner.parallelOptions;
else
    return;
end

if isempty(parallelOptions)
    return;
end
if ~isstruct(parallelOptions) || ~isscalar(parallelOptions)
    matRad_cfg.dispError('%s must be a scalar struct.',context);
end

fields = fieldnames(parallelOptions);
for i = 1:numel(fields)
    fieldName = fields{i};
    value = parallelOptions.(fieldName);
    switch fieldName
        case 'workerMemorySafetyFactor'
            options = [options {'safetyFactor',value}]; %#ok<AGROW>
        case 'memoryReserveFraction'
            options = [options {'reserveFraction',value}]; %#ok<AGROW>
        case 'minWorkerMemoryBytes'
            options = [options {'minWorkerMemoryBytes',value}]; %#ok<AGROW>
        case 'workerUpperBound'
            validateOptionalPositiveInteger(value,fieldName,matRad_cfg);
            options = [options {'workerUpperBound',value}]; %#ok<AGROW>
        otherwise
            matRad_cfg.dispError('Unknown dose parallel option "%s" in %s.', ...
                fieldName,context);
    end
end
end

function validateOptionalPositiveInteger(value,valueName,matRad_cfg)
if isempty(value)
    return;
end

if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        round(value) == value && value >= 1)
    matRad_cfg.dispError('%s must be a positive integer scalar or empty.', ...
        valueName);
end
end
