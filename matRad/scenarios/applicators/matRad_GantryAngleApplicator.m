classdef matRad_GantryAngleApplicator < matRad_ScenarioApplicator
% matRad_GantryAngleApplicator applies per-beam gantry angle uncertainty.
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
        function this = matRad_GantryAngleApplicator()
            this.name = 'gantry';
            this.componentNames = {'gantry'};
        end

        function tf = supportsComponent(~,componentName)
            tf = strncmp(componentName,'gantry.beam',numel('gantry.beam'));
        end

        function offsets = getOffsets(~,scenarioModel,scenarioId)
            offsets = scenarioModel.getGantryAngleOffset(scenarioId);
        end

        function stf = applyToStf(this,scenarioModel,scenarioId,stf)
            offsets = this.getOffsets(scenarioModel,scenarioId);
            for i = 1:min(numel(stf),numel(offsets))
                stf(i).gantryAngle = stf(i).gantryAngle + offsets(i);
            end
            stf = matRad_updateStfBeamGeometryFromAngles(stf);
        end
    end
end
