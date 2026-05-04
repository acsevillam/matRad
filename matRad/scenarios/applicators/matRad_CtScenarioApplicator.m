classdef matRad_CtScenarioApplicator < matRad_ScenarioApplicator
% matRad_CtScenarioApplicator resolves the CT scenario for a realization.
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

    methods
        function this = matRad_CtScenarioApplicator()
            this.name = 'ct';
            this.componentNames = {'ct.scenario'};
        end

        function ctScenId = getCtScenario(~,scenarioModel,scenarioId)
            ctScenId = scenarioModel.getCtScenario(scenarioId);
        end
    end
end
