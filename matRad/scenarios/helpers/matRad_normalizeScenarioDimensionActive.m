function activeDimensionNames = matRad_normalizeScenarioDimensionActive(activeDimensionNames)
% matRad_normalizeScenarioDimensionActive normalizes active uncertainty dimensions.
%
% The input can be a cell array of names, a string array, a single char
% name, or a logical mask over matRad_defaultScenarioUncertaintyDimensions().
% Returned names are ordered as ct, setup, range, gantry, couch.
%
% call
%   activeDimensionNames = matRad_normalizeScenarioDimensionActive()
%   activeDimensionNames = matRad_normalizeScenarioDimensionActive(activeDimensionNames)
%
% input
%   activeDimensionNames: active public uncertainty dimensions as names or
%                       logical mask; supported names are ct, setup, range,
%                       gantry, and couch
%
% output
%   activeDimensionNames: normalized row cell array of active uncertainty
%                       dimensions
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

matRad_cfg = MatRad_Config.instance();
validNames = matRad_defaultScenarioUncertaintyDimensions();

if nargin < 1
    activeDimensionNames = matRad_defaultScenarioActiveDimensions();
elseif isempty(activeDimensionNames)
    activeDimensionNames = {};
elseif islogical(activeDimensionNames)
    if ~isvector(activeDimensionNames) || numel(activeDimensionNames) ~= numel(validNames)
        matRad_cfg.dispError('Scenario active mask must match the number of supported uncertainty dimensions.');
    end
    activeDimensionNames = validNames(activeDimensionNames(:)');
elseif ischar(activeDimensionNames)
    activeDimensionNames = {activeDimensionNames};
elseif isstring(activeDimensionNames)
    activeDimensionNames = cellstr(activeDimensionNames(:)');
elseif iscell(activeDimensionNames)
    validCell = all(cellfun(@ischar, activeDimensionNames));
    if ~validCell
        matRad_cfg.dispError('Scenario active list must contain uncertainty dimension names.');
    end
else
    matRad_cfg.dispError('Scenario active list must be names or a logical mask.');
end

activeDimensionNames = activeDimensionNames(:)';

if numel(unique(activeDimensionNames)) ~= numel(activeDimensionNames)
    matRad_cfg.dispError('Scenario active list contains duplicate names.');
end

unknownIx = ~ismember(activeDimensionNames, validNames);
if any(unknownIx)
    matRad_cfg.dispError('Unknown scenario uncertainty dimension "%s".', activeDimensionNames{find(unknownIx, 1, 'first')});
end

[~, orderIx] = ismember(validNames, activeDimensionNames);
activeDimensionNames = validNames(orderIx > 0);

end
