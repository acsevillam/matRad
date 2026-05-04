function [scenarioValues,ctScenIds,scenProb,scenWeight,scenForProb,linearMask,scenMask] = matRad_filterZeroProbabilityScenarios(scenarioValues,ctScenIds,scenProb,scenWeight,scenForProb,linearMask,scenMask)
% matRad_filterZeroProbabilityScenarios removes scenarios with zero weight.
%
% Probability-weighted scenario models should pass scenProb as scenWeight
% before normalization. Sampled models can pass their sampling weights.
% scenProb is preserved as the unnormalized probability/density value, while
% scenWeight is filtered and normalized to sum to one.
%
% call
%   [scenarioValues,ctScenIds,scenProb,scenWeight,scenForProb,linearMask,scenMask] = matRad_filterZeroProbabilityScenarios(scenarioValues,ctScenIds,scenProb,scenWeight,scenForProb,linearMask,scenMask)
%
% input
%   scenarioValues:      scenario realization values, one row per scenario
%   ctScenIds:           CT scenario id for each realization row
%   scenProb:            unnormalized probability/density values
%   scenWeight:          scenario weights used for optimization/statistics
%   scenForProb:         legacy probability realization matrix
%   linearMask:          scenario mask subscripts, one row per scenario
%   scenMask:            logical DIJ scenario mask
%
% output
%   scenarioValues:      filtered scenario realization values
%   ctScenIds:           filtered CT scenario ids
%   scenProb:            filtered unnormalized probability/density values
%   scenWeight:          filtered and normalized scenario weights
%   scenForProb:         filtered legacy probability realization matrix
%   linearMask:          filtered scenario mask subscripts
%   scenMask:            updated logical DIJ scenario mask
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

numScenarios = size(scenarioValues,1);
ctScenIds = ctScenIds(:);
scenProb = scenProb(:);
scenWeight = scenWeight(:);

validRows = numel(ctScenIds) == numScenarios && ...
    numel(scenProb) == numScenarios && ...
    numel(scenWeight) == numScenarios && ...
    size(scenForProb,1) == numScenarios && ...
    size(linearMask,1) == numScenarios;

if ~validRows
    matRad_cfg.dispError('Scenario realization arrays must have consistent row counts.');
end

validProbabilityValues = all(isfinite(scenProb)) && all(scenProb >= 0) && ...
    all(isfinite(scenWeight)) && all(scenWeight >= 0);

if ~validProbabilityValues
    matRad_cfg.dispError('Scenario probabilities and weights must be finite non-negative values.');
end

keepIx = scenWeight > 0;

if ~any(keepIx)
    matRad_cfg.dispError('Scenario probabilities must contain at least one positive active scenario.');
end

if ~all(keepIx)
    scenarioValues = scenarioValues(keepIx,:);
    ctScenIds = ctScenIds(keepIx);
    scenProb = scenProb(keepIx);
    scenWeight = scenWeight(keepIx);
    scenForProb = scenForProb(keepIx,:);
    linearMask = linearMask(keepIx,:);

    scenMask = false(size(scenMask));
    if ~isempty(linearMask)
        scenMaskSize = size(scenMask);
        if numel(scenMaskSize) < size(linearMask,2)
            scenMaskSize(end+1:size(linearMask,2)) = 1;
        end
        maskSubscripts = mat2cell(linearMask, size(linearMask,1), ones(1,size(linearMask,2)));
        maskIx = sub2ind(scenMaskSize,maskSubscripts{:});
        scenMask(maskIx) = true;
    end
end

scenWeight = scenWeight ./ sum(scenWeight);

end
