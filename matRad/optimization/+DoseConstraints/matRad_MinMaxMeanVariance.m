classdef matRad_MinMaxMeanVariance < DoseConstraints.matRad_DoseConstraint
% matRad_MinMaxMeanVariance Implements a PROB2 mean variance constraint
%   See matRad_DoseConstraint for interface description
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
        name = 'mean variance constraint';
        parameterNames = {'sigma^2_{min}','sigma^2_{max}'};
        parameterTypes = {'numeric','numeric'};
    end

    properties
        parameters = {0,Inf};
    end

    methods
        function this = matRad_MinMaxMeanVariance(minMeanVariance,maxMeanVariance)
            if nargin == 1 && isstruct(minMeanVariance)
                inputStruct = minMeanVariance;
                initFromStruct = true;
            else
                inputStruct = [];
                initFromStruct = false;
            end

            this@DoseConstraints.matRad_DoseConstraint(inputStruct);

            if ~initFromStruct
                if nargin >= 2 && isscalar(maxMeanVariance)
                    this.parameters{2} = maxMeanVariance;
                end
                if nargin >= 1 && isscalar(minMeanVariance)
                    this.parameters{1} = minMeanVariance;
                end
            end
        end

        function cu = upperBounds(this,~)
            cu = this.parameters{2};
        end

        function cl = lowerBounds(this,~)
            cl = this.parameters{1};
        end

        function cDose = computeDoseConstraintFunction(~,~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError(['Mean Variance constraint is only ', ...
                'supported for PROB2 fluence optimization!']);
            cDose = [];
        end

        function cDoseJacob = computeDoseConstraintJacobian(~,~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError(['Mean Variance constraint has no ', ...
                'nominal dose-jacobian implementation!']);
            cDoseJacob = [];
        end

        function c = computeProb2ConstraintFunction(~,stats)
            c = stats.meanVariance;
        end

        function jacob = computeProb2ConstraintJacobian(~,stats)
            jacob = stats.gradMeanVariance(:);
        end
    end

    methods (Static)
        function rob = availableRobustness()
            rob = {'PROB2'};
        end
    end
end
