classdef matRad_SetupShiftApplicator < matRad_ScenarioApplicator
    % matRad_SetupShiftApplicator applies setup uncertainty realizations.
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

        function this = matRad_SetupShiftApplicator()
            this.name = 'setup';
            this.componentNames = {'setup.x', 'setup.y', 'setup.z'};
        end

        function shift = getShift(~, scenarioModel, scenarioId)
            shift = scenarioModel.getSetupShift(scenarioId);
        end

        function stf = applyToStf(this, scenarioModel, scenarioId, stf)
            shift = this.getShift(scenarioModel, scenarioId);
            for i = 1:numel(stf)
                stf(i).isoCenter = stf(i).isoCenter + shift;
            end
        end

        function margin = getMargin(~, scenarioModel, pbMargin)
            if nargin < 3 || isempty(pbMargin)
                pbMargin = 0;
            end
            scenarioIds = scenarioModel.scenarioIds();
            shifts = zeros(numel(scenarioIds), 3);
            for i = 1:numel(scenarioIds)
                shifts(i, :) = scenarioModel.getSetupShift(scenarioIds(i));
            end
            maxShift = max(abs(shifts), [], 1);
            margin = maxShift + pbMargin;
        end

    end
end
