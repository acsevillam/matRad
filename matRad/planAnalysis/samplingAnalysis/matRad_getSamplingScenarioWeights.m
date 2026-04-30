function scenWeights = matRad_getSamplingScenarioWeights(pln,numScenarios,ctScenProb,inputScenWeights)
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
%                       [ctScenIndex probability]
%   inputScenWeights:   (optional) resolved scenario weights to use directly
%
% output
%   scenWeights:        normalized scenario weights for all sampled scenarios
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

if nargin < 3
    ctScenProb = [];
end

if nargin < 4
    inputScenWeights = [];
end

if ~isempty(inputScenWeights)
    scenWeights = inputScenWeights(:);
    validateScenarioWeightCount(scenWeights,numScenarios);
    scenWeights = normalizeScenarioWeights(scenWeights);
    return;
end

if hasFieldOrProp(pln.multScen,'scenWeight') && ~isempty(pln.multScen.scenWeight)
    scenWeights = pln.multScen.scenWeight(:);
elseif hasFieldOrProp(pln.multScen,'scenProb') && ~isempty(pln.multScen.scenProb)
    scenWeights = pln.multScen.scenProb(:);
else
    matRad_cfg.dispError('Scenario model must provide scenWeight or scenProb.');
end

validateScenarioWeightCount(scenWeights,numScenarios);

if isempty(ctScenProb)
    scenWeights = normalizeScenarioWeights(scenWeights);
    return;
end

validateCtScenProb(ctScenProb);

if ~hasFieldOrProp(pln.multScen,'linearMask') || isempty(pln.multScen.linearMask)
    matRad_cfg.dispError('ctScenProb override requires multScen.linearMask.');
end

if ~hasFieldOrProp(pln.multScen,'ctScenProb') || isempty(pln.multScen.ctScenProb)
    matRad_cfg.dispError('ctScenProb override requires multScen.ctScenProb.');
end

ctScen = pln.multScen.linearMask(:,1);
if numel(ctScen) ~= numScenarios
    matRad_cfg.dispError('Number of CT scenario indices does not match number of sampled scenarios.');
end

oldCtProb = getCtScenarioProbabilities(pln.multScen.ctScenProb,ctScen);
newCtProb = getCtScenarioProbabilities(ctScenProb,ctScen);

invalidOriginalProb = oldCtProb <= 0 & scenWeights > 0;
if any(invalidOriginalProb)
    matRad_cfg.dispError('Cannot reweight scenarios with positive weight and zero original CT scenario probability.');
end

oldCtProb(oldCtProb <= 0) = 1;
scenWeights = scenWeights ./ oldCtProb .* newCtProb;
scenWeights = normalizeScenarioWeights(scenWeights);

end

function validateScenarioWeightCount(scenWeights,numScenarios)
if numel(scenWeights) ~= numScenarios
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Number of scenario weights does not match number of sampled scenarios.');
end
end

function validateCtScenProb(ctScenProb)
valid = isnumeric(ctScenProb) && ismatrix(ctScenProb) && size(ctScenProb,2) == 2 && ...
    all(isfinite(ctScenProb(:))) && all(ctScenProb(:,2) >= 0) && ...
    all(round(ctScenProb(:,1)) == ctScenProb(:,1)) && all(ctScenProb(:,1) > 0);
if ~valid
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('ctScenProb must be a two-column matrix [ctScenIndex probability].');
end

if sum(ctScenProb(:,2)) <= 0
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('ctScenProb must contain at least one positive probability.');
end
end

function ctProb = getCtScenarioProbabilities(ctScenProb,ctScen)
matRad_cfg = MatRad_Config.instance();
ctProb = zeros(size(ctScen));

for i = 1:numel(ctScen)
    ctProbIx = find(ctScenProb(:,1) == ctScen(i),1,'first');
    if isempty(ctProbIx)
        matRad_cfg.dispError('Missing probability for CT scenario %d.',ctScen(i));
    end
    ctProb(i) = ctScenProb(ctProbIx,2);
end
end

function scenWeights = normalizeScenarioWeights(scenWeights)
matRad_cfg = MatRad_Config.instance();

scenWeights = scenWeights(:);
scenWeights(~isfinite(scenWeights) | scenWeights < 0) = 0;
if sum(scenWeights) <= 0
    matRad_cfg.dispError('Scenario weights must contain at least one positive finite value.');
end
scenWeights = scenWeights./sum(scenWeights);
end

function tf = hasFieldOrProp(value,fieldName)
tf = (isobject(value) && isprop(value,fieldName)) || ...
    (isstruct(value) && isfield(value,fieldName));
end
