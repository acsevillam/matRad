function scenProb = matRad_computeScenarioProbabilities(components, scenarioValues, ctScenProb, scenarioCtScenIds)
% matRad_computeScenarioProbabilities evaluates scenario probabilities.
%
% Scenario probabilities are evaluated as normal-density values at the
% scenario realization positions and multiplied by CT scenario probabilities.
% Inactive components are ignored, which keeps inactive or zero-uncertainty
% dimensions from making the covariance matrix singular.
%
% call
%   scenProb = matRad_computeScenarioProbabilities(components, scenarioValues, ctScenProb, scenarioCtScenIds)
%
% input
%   components:          scenario component struct array with active flags
%                       and uncertainty scales
%   scenarioValues:      matrix with one scenario realization per row and
%                       one component value per column
%   ctScenProb:          Nx2 matrix; first column is CT scenario id and
%                       second column is CT scenario probability
%   scenarioCtScenIds:   CT scenario id for each scenarioValues row
%
% output
%   scenProb:            unnormalized probability/density value for each
%                       scenario realization row
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

if isempty(scenarioValues)
    scenProb = zeros(0, 1);
    return
end

helper_validateInputs(components, scenarioValues, ctScenProb, scenarioCtScenIds);

activeIx = logical([components.active]);
if any(activeIx)
    scales = [components(activeIx).scale];
    activeValues = scenarioValues(:, activeIx);
    sigma = diag(scales.^2);
    numComponents = size(sigma, 1);
    cholSigma = chol(sigma);
    realizationProb = (2 * pi)^(-numComponents / 2) .* ...
        exp(-0.5 .* sum((activeValues / cholSigma).^2, 2)) ./ prod(diag(cholSigma));
else
    realizationProb = ones(size(scenarioValues, 1), 1);
end

phaseProb = zeros(numel(scenarioCtScenIds), 1);
for i = 1:numel(scenarioCtScenIds)
    ctScenIx = find(ctScenProb(:, 1) == scenarioCtScenIds(i), 1, 'first');
    helper_errorIf(isempty(ctScenIx), 'Could not find CT scenario %d in ctScenProb.', scenarioCtScenIds(i));
    phaseProb(i) = ctScenProb(ctScenIx, 2);
end

scenProb = phaseProb(:) .* realizationProb(:);

end

function helper_validateInputs(components, scenarioValues, ctScenProb, scenarioCtScenIds)

helper_validateComponents(components);
helper_validateScenarioValues(scenarioValues, components, scenarioCtScenIds);
helper_validateCtScenProb(ctScenProb);
helper_validateCtScenIds(scenarioCtScenIds);

end

function helper_validateComponents(components)

requiredFields = {'name', 'applicator', 'unit', 'active', 'scale'};
validComponents = isstruct(components) && all(isfield(components, requiredFields));
helper_errorIf(~validComponents, ...
               'Scenario components must be a struct array with name, applicator, unit, active, and scale fields.');

for i = 1:numel(components)
    component = components(i);
    helper_errorIf(~helper_isValidComponent(component), ...
                   'Scenario component metadata is invalid.');

    helper_errorIf(all([logical(component.active), component.scale <= 0]), ...
                   'Active scenario component "%s" requires a positive uncertainty scale.', component.name);
end

end

function tf = helper_isValidComponent(component)

tf = all([ischar(component.name), ...
          ~isempty(component.name), ...
          isscalar(component.active), ...
          islogical(component.active) || isnumeric(component.active), ...
          isnumeric(component.scale), ...
          isscalar(component.scale), ...
          isfinite(component.scale), ...
          component.scale >= 0]);

end

function helper_validateScenarioValues(scenarioValues, components, scenarioCtScenIds)

helper_errorIf(~(isnumeric(scenarioValues) && all(isfinite(scenarioValues(:)))), ...
               'Scenario values must be finite numeric values.');
helper_errorIf(size(scenarioValues, 2) ~= numel(components), ...
               'Scenario values must contain one column per scenario component.');
helper_errorIf(numel(scenarioCtScenIds) ~= size(scenarioValues, 1), ...
               'Scenario CT scenario ids must contain one entry per scenario realization.');

end

function helper_validateCtScenProb(ctScenProb)

validShape = isnumeric(ctScenProb) && ismatrix(ctScenProb) && ...
    ~isempty(ctScenProb) && size(ctScenProb, 2) == 2;

validValues = false;
if validShape
    validValues = all(isfinite(ctScenProb(:))) && ...
        all(round(ctScenProb(:, 1)) == ctScenProb(:, 1)) && ...
        all(ctScenProb(:, 1) >= 1) && all(ctScenProb(:, 2) > 0);
end

helper_errorIf(~(validShape && validValues), ...
               'ctScenProb must contain positive integer CT scenario ids and positive probabilities.');

end

function helper_validateCtScenIds(scenarioCtScenIds)

validCtScenIds = false;
if isnumeric(scenarioCtScenIds)
    validCtScenIds = all(isfinite(scenarioCtScenIds(:))) && ...
        all(round(scenarioCtScenIds(:)) == scenarioCtScenIds(:)) && ...
        all(scenarioCtScenIds(:) >= 1);
end

helper_errorIf(~validCtScenIds, ...
               'Scenario CT scenario ids must be positive integer values.');

end

function helper_errorIf(condition, message, varargin)

if condition
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError(message, varargin{:});
end

end
