classdef matRad_MeanVariance < DoseObjectives.matRad_DoseObjective
% matRad_MeanVariance Implements a PROB2 mean variance objective
%   See matRad_DoseObjective for interface description
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

    properties (Constant)
        name = 'Mean Variance';
        parameterNames = {};
        parameterTypes = {};
    end

    properties
        parameters = {};
        penalty = 1;
    end

    methods
        function obj = matRad_MeanVariance(penalty)
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                inputStruct = [];
                initFromStruct = false;
            end

            obj@DoseObjectives.matRad_DoseObjective(inputStruct);

            if ~initFromStruct && nargin >= 1 && isscalar(penalty)
                obj.penalty = penalty;
            end
        end

        function fDose = computeDoseObjectiveFunction(~,~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError(['Mean Variance objective is only ', ...
                'supported for PROB2 fluence optimization!']);
            fDose = [];
        end

        function fDoseGrad = computeDoseObjectiveGradient(~,~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError(['Mean Variance objective has no ', ...
                'nominal dose-gradient implementation!']);
            fDoseGrad = [];
        end

        function f = computeProb2ObjectiveFunction(~,stats)
            f = stats.meanVariance;
        end

        function fGrad = computeProb2ObjectiveGradient(~,stats)
            fGrad = stats.gradMeanVariance;
        end
    end

    methods (Static)
        function rob = availableRobustness()
            rob = {'PROB2'};
        end
    end
end
