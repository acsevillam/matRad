classdef matRad_SquaredBertoluzzaDeviation < DoseObjectives.matRad_DoseObjective
    % matRad_SquaredBertoluzzaDeviation Implements a squared Bertoluzza objective
    %   See matRad_DoseObjective for interface description.
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
        name = 'Squared Bertoluzza Deviation'
        parameterNames = {'d^{ref}'}
        parameterTypes = {'dose'}
    end

    properties
        parameters = {60}
        penalty = 1
    end

    methods

        function obj = matRad_SquaredBertoluzzaDeviation(penalty, dRef)
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                inputStruct = [];
                initFromStruct = false;
            end

            obj@DoseObjectives.matRad_DoseObjective(inputStruct);

            if ~initFromStruct
                if nargin >= 2 && isscalar(dRef)
                    obj.parameters{1} = dRef;
                end

                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end

        function fDose = computeDoseObjectiveFunction(obj, varargin)
            if numel(varargin) == 4
                fDose = obj.computeBertoluzzaObjective(varargin{:});
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError(['Squared Bertoluzza objective is only ', ...
                                      'supported for INTERVAL2/INTERVAL3 ', ...
                                      'fluence optimization!']);
            end
        end

        function fDoseGrad = computeDoseObjectiveGradient(~, ~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError(['Squared Bertoluzza objective has no ', ...
                                  'nominal dose-gradient implementation!']);
            fDoseGrad = [];
        end

        function fWGrad = computeFluenceObjectiveGradient(obj, w, ix, theta, dijInterval)
            [centerRows, radiusMatrix] = obj.getIntervalMatrices(w, ix, dijInterval);
            centerDose = centerRows * w;
            deviation = centerDose - obj.parameters{1};

            centerGradient = 2 * (centerRows' * deviation);
            radiusGradient = (radiusMatrix + radiusMatrix') * w;
            centerNormGradient = 2 * (centerRows' * centerDose);

            fWGrad = obj.penalty / numel(ix) * ...
                (centerGradient + theta * (radiusGradient - centerNormGradient));
        end

    end

    methods (Access = private)

        function fFluence = computeBertoluzzaObjective(obj, w, ix, theta, dijInterval)
            [centerRows, radiusMatrix] = obj.getIntervalMatrices(w, ix, dijInterval);
            centerDose = centerRows * w;
            deviation = centerDose - obj.parameters{1};

            radiusTerm = w' * radiusMatrix * w;
            centerTerm = centerDose' * centerDose;

            fFluence = obj.penalty / numel(ix) * ...
                (deviation' * deviation + theta * (radiusTerm - centerTerm));
        end

        function [centerRows, radiusMatrix] = getIntervalMatrices(~, w, ix, dijInterval)
            matRad_cfg = MatRad_Config.instance();

            if ~isfield(dijInterval, 'center') || ~isfield(dijInterval, 'radius')
                matRad_cfg.dispError('dij_interval must contain center and radius fields!');
            end

            centerMatrix = dijInterval.center;
            radiusMatrix = dijInterval.radius;
            ix = ix(:);

            if isempty(ix)
                matRad_cfg.dispError('Squared Bertoluzza objective requires non-empty voxel indices!');
            end

            if size(centerMatrix, 2) ~= numel(w) || ...
                    size(radiusMatrix, 1) ~= numel(w) || ...
                    size(radiusMatrix, 2) ~= numel(w)
                matRad_cfg.dispError('dij_interval dimensions are inconsistent with the fluence vector!');
            end

            if any(ix < 1) || any(ix > size(centerMatrix, 1))
                matRad_cfg.dispError('Voxel indices exceed dij_interval.center dimensions!');
            end

            centerRows = centerMatrix(ix, :);
        end

    end

    methods (Static)

        function rob = availableRobustness()
            rob = {'INTERVAL2', 'INTERVAL3'};
        end

    end
end
