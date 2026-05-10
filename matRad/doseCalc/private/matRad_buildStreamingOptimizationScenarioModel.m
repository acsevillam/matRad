function scenarioModel = matRad_buildStreamingOptimizationScenarioModel(refScen)
% matRad_buildStreamingOptimizationScenarioModel creates a nominal model for streaming contexts
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

components = matRad_createScenarioComponents([0 0 0],0,0,{'ct'});
scenarioValues = zeros(1,numel(components));
ctScenIds = refScen;
scenProb = 1;
scenWeight = 1;
scenForProb = [ctScenIds scenarioValues];
linearMask = [1 1 1];
scenMask = true(1,1,1);

scenarioModel = matRad_NominalScenario();
scenarioModel.setScenarioRealizations(components,scenarioValues,ctScenIds, ...
    scenProb,scenWeight,scenForProb,linearMask,scenMask);
end
