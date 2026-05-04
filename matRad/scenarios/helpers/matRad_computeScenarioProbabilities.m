function scenProb = matRad_computeScenarioProbabilities(components,scenarioValues,ctScenProb,scenarioCtScenIds)
% matRad_computeScenarioProbabilities evaluates scenario probabilities.
%
% Components belonging to inactive dimensions/applicators are ignored. This
% keeps zero-uncertainty components exact instead of replacing them with eps.
% Continuous realization probabilities are evaluated as normal-density
% values at the scenario grid/sample positions and multiplied by CT phase
% probabilities. They are not normalized scenario weights; callers normalize
% scenWeight separately when needed.
%
% call
%   scenProb = matRad_computeScenarioProbabilities(components,scenarioValues,ctScenProb,scenarioCtScenIds)
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
    scenProb = [];
    return;
end

activeIx = [components.active];
if any(activeIx)
    scales = [components(activeIx).scale];
    activeValues = scenarioValues(:,activeIx);
    Sigma = diag(scales.^2);
    d = size(Sigma,1);
    cs = chol(Sigma);
    realizationProb = (2*pi)^(-d/2) * ...
        exp(-0.5*sum((activeValues/cs).^2,2)) / prod(diag(cs));
else
    realizationProb = ones(size(scenarioValues,1),1);
end

phaseProb = zeros(numel(scenarioCtScenIds),1);
for i = 1:numel(scenarioCtScenIds)
    ctScenIx = find(ctScenProb(:,1) == scenarioCtScenIds(i),1,'first');
    if isempty(ctScenIx)
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Could not find CT scenario %d in ctScenProb.',scenarioCtScenIds(i));
    end
    phaseProb(i) = ctScenProb(ctScenIx,2);
end

scenProb = phaseProb(:) .* realizationProb(:);

end
