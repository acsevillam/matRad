function scenarioModel = matRad_buildNominalScenarioModel(ctScenId)
% matRad_buildNominalScenarioModel builds one CT-only nominal scenario model
%
% call
%   scenarioModel = ScenarioBatch.Scenarios.matRad_buildNominalScenarioModel(ctScenId)
%
% input
%   ctScenId:       CT scenario id represented by the nominal realization
%
% output
%   scenarioModel:  compact one-realization matRad_ScenarioModel
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

matRadCfg = MatRad_Config.instance();
if nargin < 1 || isempty(ctScenId)
    ctScenId = 1;
end
if ~isnumeric(ctScenId) || ~isscalar(ctScenId) || ~isfinite(ctScenId) || ...
        ctScenId < 1 || fix(ctScenId) ~= ctScenId
    matRadCfg.dispError('Nominal scenario CT id must be a positive integer.');
end

components = matRad_createScenarioComponents([0 0 0], 0, 0, {'ct'});
scenarioValues = zeros(1, numel(components));
ctScenIds = double(ctScenId);
scenProb = 1;
scenWeight = 1;
scenForProb = [ctScenIds scenarioValues];
storageSubscripts = [1 1];
storageSize = [1 1];

scenarioModel = matRad_NominalScenario();
scenarioModel.setScenarioRealizations(components, scenarioValues, ctScenIds, ...
                                      scenProb, scenWeight, scenForProb, ...
                                      storageSubscripts, storageSize, 'compact-realization');
end
