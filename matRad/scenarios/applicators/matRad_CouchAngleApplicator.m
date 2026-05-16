classdef matRad_CouchAngleApplicator < matRad_ScenarioApplicator
    % matRad_CouchAngleApplicator applies per-beam couch angle uncertainty.
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

        function this = matRad_CouchAngleApplicator()
            this.name = 'couch';
            this.componentNames = {'couch'};
        end

        function tf = supportsComponent(~, componentName)
            tf = strncmp(componentName, 'couch.beam', numel('couch.beam'));
        end

        function offsets = getOffsets(~, scenarioModel, scenarioId)
            offsets = scenarioModel.getCouchAngleOffset(scenarioId);
        end

        function stf = applyToStf(this, scenarioModel, scenarioId, stf)
            offsets = this.getOffsets(scenarioModel, scenarioId);
            for i = 1:min(numel(stf), numel(offsets))
                stf(i).couchAngle = stf(i).couchAngle + offsets(i);
            end
            stf = matRad_updateStfBeamGeometryFromAngles(stf);
        end

    end
end
