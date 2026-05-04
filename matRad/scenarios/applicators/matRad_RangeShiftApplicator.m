classdef matRad_RangeShiftApplicator < matRad_ScenarioApplicator
% matRad_RangeShiftApplicator applies range uncertainty realizations.
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
        function this = matRad_RangeShiftApplicator()
            this.name = 'range';
            this.componentNames = {'range.absolute','range.relative'};
        end

        function [absRangeShift,relRangeShift] = getShift(~,scenarioModel,scenarioId)
            rangeShift = scenarioModel.getRangeShift(scenarioId);
            absRangeShift = rangeShift(1);
            relRangeShift = rangeShift(2);
        end

        function radDepths = applyToRadDepths(this,scenarioModel,scenarioId,radDepths)
            [absRangeShift,relRangeShift] = this.getShift(scenarioModel,scenarioId);
            radDepths = (1+relRangeShift) .* radDepths + absRangeShift;
            if absRangeShift < 0
                radDepths(radDepths < 0) = 0;
            end
        end

        function radDepths = applyToCumulativeDepths(this,scenarioModel,scenarioId,radDepths,rho)
            [absRangeShift,relRangeShift] = this.getShift(scenarioModel,scenarioId);
            if relRangeShift ~= 0 || absRangeShift ~= 0
                radDepths = radDepths + rho .* relRangeShift + absRangeShift;
                radDepths(radDepths < 0) = 0;
            end
        end

        function margin = getMargin(~,scenarioModel)
            scenarioIds = scenarioModel.scenarioIds();
            rangeShifts = zeros(numel(scenarioIds),2);
            for i = 1:numel(scenarioIds)
                rangeShifts(i,:) = scenarioModel.getRangeShift(scenarioIds(i));
            end
            margin = max(abs(rangeShifts(:,1))) + max(abs(rangeShifts(:,2)));
        end
    end
end
