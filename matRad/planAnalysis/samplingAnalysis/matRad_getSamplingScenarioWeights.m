function scenWeights = matRad_getSamplingScenarioWeights(pln, numScenarios, ctScenProb, inputScenWeights)
% matRad_getSamplingScenarioWeights resolves sampling scenario weights
%
% call
%   scenWeights = matRad_getSamplingScenarioWeights(pln,numScenarios)
%   scenWeights = matRad_getSamplingScenarioWeights(pln,numScenarios,ctScenProb)
%   scenWeights = matRad_getSamplingScenarioWeights(pln,numScenarios,ctScenProb,inputScenWeights)
%
% input
%   pln:                matRad pln struct with multScen scenario model
%   numScenarios:       number of sampled scenarios
%   ctScenProb:         (optional) CT scenario probability override as
%                       [ctScenId probability]
%   inputScenWeights:   (optional) resolved scenario weights to use directly
%
% output
%   scenWeights:        normalized scenario weights for all sampled scenarios
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

if nargin < 3
    ctScenProb = [];
end

if nargin < 4
    inputScenWeights = [];
end

if ~isfield(pln, 'multScen') || isempty(pln.multScen)
    matRad_cfg.dispError('Sampling scenario weights require pln.multScen.');
end

if ~isempty(inputScenWeights)
    scenWeights = inputScenWeights(:);
    matRad_validateScenarioWeightCount(scenWeights, numScenarios);
    scenWeights = matRad_normalizeScenarioWeights(scenWeights);
    return
end

if matRad_hasFieldOrProp(pln.multScen, 'scenWeight') && ~isempty(pln.multScen.scenWeight)
    scenWeights = pln.multScen.scenWeight(:);
elseif matRad_hasFieldOrProp(pln.multScen, 'scenProb') && ~isempty(pln.multScen.scenProb)
    scenWeights = pln.multScen.scenProb(:);
else
    matRad_cfg.dispError('Scenario model must provide scenWeight or scenProb.');
end

matRad_validateScenarioWeightCount(scenWeights, numScenarios);

if isempty(ctScenProb)
    scenWeights = matRad_normalizeScenarioWeights(scenWeights);
    return
end

matRad_validateCtScenProb(ctScenProb);

if ~isa(pln.multScen, 'matRad_ScenarioModel')
    matRad_cfg.dispError('ctScenProb override requires a matRad_ScenarioModel instance.');
end

if ~matRad_hasFieldOrProp(pln.multScen, 'ctScenProb') || isempty(pln.multScen.ctScenProb)
    matRad_cfg.dispError('ctScenProb override requires multScen.ctScenProb.');
end

scenarioIds = pln.multScen.scenarioIds();
ctScenIds = arrayfun(@(id) pln.multScen.getCtScenario(id), scenarioIds);
if numel(ctScenIds) ~= numScenarios
    matRad_cfg.dispError('Number of CT scenario ids does not match number of sampled scenarios.');
end

oldCtProb = matRad_getCtScenarioProbabilities(pln.multScen.ctScenProb, ctScenIds);
newCtProb = matRad_getCtScenarioProbabilities(ctScenProb, ctScenIds);

invalidOriginalProb = oldCtProb <= 0 & scenWeights > 0;
if any(invalidOriginalProb)
    matRad_cfg.dispError('Cannot reweight scenarios with positive weight and zero original CT scenario probability.');
end

oldCtProb(oldCtProb <= 0) = 1;
scenWeights = scenWeights ./ oldCtProb .* newCtProb;
scenWeights = matRad_normalizeScenarioWeights(scenWeights);

end

function matRad_validateScenarioWeightCount(scenWeights, numScenarios)
if numel(scenWeights) ~= numScenarios
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Number of scenario weights does not match number of sampled scenarios.');
end
end

function matRad_validateCtScenProb(ctScenProb)
valid = isnumeric(ctScenProb) && ismatrix(ctScenProb) && size(ctScenProb, 2) == 2 && ...
    all(isfinite(ctScenProb(:))) && all(ctScenProb(:, 2) >= 0) && ...
    all(round(ctScenProb(:, 1)) == ctScenProb(:, 1)) && all(ctScenProb(:, 1) > 0);
if ~valid
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('ctScenProb must be a two-column matrix [ctScenId probability].');
end

if sum(ctScenProb(:, 2)) <= 0
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('ctScenProb must contain at least one positive probability.');
end
end

function ctProb = matRad_getCtScenarioProbabilities(ctScenProb, ctScenIds)
matRad_cfg = MatRad_Config.instance();
ctProb = zeros(size(ctScenIds));

for i = 1:numel(ctScenIds)
    ctScenIx = find(ctScenProb(:, 1) == ctScenIds(i), 1, 'first');
    if isempty(ctScenIx)
        matRad_cfg.dispError('Missing probability for CT scenario %d.', ctScenIds(i));
    end
    ctProb(i) = ctScenProb(ctScenIx, 2);
end
end

function scenWeights = matRad_normalizeScenarioWeights(scenWeights)
matRad_cfg = MatRad_Config.instance();

scenWeights = scenWeights(:);
scenWeights(~isfinite(scenWeights) | scenWeights < 0) = 0;
if sum(scenWeights) <= 0
    matRad_cfg.dispError('Scenario weights must contain at least one positive finite value.');
end
scenWeights = scenWeights ./ sum(scenWeights);
end

function tf = matRad_hasFieldOrProp(value, fieldName)
tf = (isobject(value) && isprop(value, fieldName)) || ...
    (isstruct(value) && isfield(value, fieldName));
end
