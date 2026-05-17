classdef matRad_MinMaxMeanVariance < DoseConstraints.matRad_DoseConstraint
    % matRad_MinMaxMeanVariance Implements a PROB2 mean variance constraint
    %   See matRad_DoseConstraint for interface description.
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
        name = 'Mean Variance'
        parameterNames = {'sigma^2_{min}', 'sigma^2_{max}'}
        parameterTypes = {'numeric', 'numeric'}
    end

    properties
        parameters = {0, Inf}
    end

    methods

        function obj = matRad_MinMaxMeanVariance(minMeanVariance, maxMeanVariance)
            if nargin == 1 && isstruct(minMeanVariance)
                inputStruct = minMeanVariance;
                initFromStruct = true;
            else
                inputStruct = [];
                initFromStruct = false;
            end

            obj@DoseConstraints.matRad_DoseConstraint(inputStruct);

            if ~initFromStruct
                if nargin >= 2 && isscalar(maxMeanVariance)
                    obj.parameters{2} = maxMeanVariance;
                end
                if nargin >= 1 && isscalar(minMeanVariance)
                    obj.parameters{1} = minMeanVariance;
                end
            end
        end

        function cu = upperBounds(obj, ~)
            cu = obj.parameters{2};
        end

        function cl = lowerBounds(obj, ~)
            cl = obj.parameters{1};
        end

        function cDose = computeDoseConstraintFunction(~, ~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError(['Mean Variance constraint is only ', ...
                                  'supported for PROB2 fluence optimization!']);
            cDose = [];
        end

        function cDoseJacob = computeDoseConstraintJacobian(~, ~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError(['Mean Variance constraint has no ', ...
                                  'nominal dose-jacobian implementation!']);
            cDoseJacob = [];
        end

        function c = computeProb2ConstraintFunction(~, stats)
            if isempty(stats.meanVariance)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('PROB2 mean variance requires an Omega matrix for the selected structure!');
            end
            c = stats.meanVariance;
        end

        function jacob = computeProb2ConstraintJacobian(~, stats)
            if isempty(stats.gradMeanVariance)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('PROB2 mean variance requires an Omega matrix for the selected structure!');
            end
            jacob = stats.gradMeanVariance(:);
        end

    end

    methods (Static)

        function rob = availableRobustness()
            rob = {'PROB2'};
        end

    end

end
