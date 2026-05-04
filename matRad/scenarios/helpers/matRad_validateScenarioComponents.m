function matRad_validateScenarioComponents(components)
% matRad_validateScenarioComponents validates scenario realization components.
%
% Active continuous components must have a strictly positive scale. Inactive
% components stay present for compatibility, but are ignored by scenario
% generation and probability evaluation through their parent applicator.
% Required component fields are name, applicator, unit, active, and scale.
%
% call
%   matRad_validateScenarioComponents(components)
%
% input
%   components:          scenario component struct array to validate
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

requiredFields = {'name','applicator','unit','active','scale'};
validFields = isstruct(components) && all(isfield(components,requiredFields));
if ~validFields
    matRad_cfg.dispError('Scenario components must be a struct array with name, applicator, unit, active, and scale fields.');
end

names = {components.name};
if numel(unique(names)) ~= numel(names)
    matRad_cfg.dispError('Scenario component names must be unique.');
end

for i = 1:numel(components)
    component = components(i);
    validName = ischar(component.name) && ~isempty(component.name);
    validActive = isscalar(component.active) && ...
        (islogical(component.active) || isnumeric(component.active));
    validScale = isnumeric(component.scale) && isscalar(component.scale) && ...
        isfinite(component.scale) && component.scale >= 0;

    if ~validName || ~validActive || ~validScale
        matRad_cfg.dispError('Scenario component "%s" has invalid metadata.',component.name);
    end

    if logical(component.active) && component.scale <= 0
        matRad_cfg.dispError('Active scenario component "%s" requires a positive uncertainty scale.',component.name);
    end
end

end
