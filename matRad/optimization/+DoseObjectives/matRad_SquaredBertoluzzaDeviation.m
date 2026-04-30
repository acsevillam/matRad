classdef matRad_SquaredBertoluzzaDeviation < DoseObjectives.matRad_DoseObjective
% matRad_SquaredBertoluzzaDeviation Implements a squared Bertoluzza objective
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
        name = 'Squared Bertoluzza Deviation';
        parameterNames = {'d^{ref}'};
        parameterTypes = {'dose'};
    end

    properties
        parameters = {60};
        penalty = 1;
    end

    methods
        function obj = matRad_SquaredBertoluzzaDeviation(penalty,dRef)

            %If we have a struct in first argument
            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end

            %Call Superclass Constructor (for struct initialization)
            obj@DoseObjectives.matRad_DoseObjective(inputStruct);

            %now handle initialization from other parameters
            if ~initFromStruct
                if nargin >= 2 && isscalar(dRef)
                    obj.parameters{1} = dRef;
                end

                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end

        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,varargin)
            if numel(varargin) == 4
                fDose = computeBertoluzzaObjective(obj,varargin{:});
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError(['Squared Bertoluzza objective is only ', ...
                    'supported for INTERVAL2/INTERVAL3 fluence optimization!']);
            end
        end

        %% Calculates the Objective Function gradient
        function fDoseGrad = computeDoseObjectiveGradient(~,~)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError(['Squared Bertoluzza objective has no ', ...
                'nominal dose-gradient implementation!']);
            fDoseGrad = [];
        end

        function fWGrad = computeFluenceObjectiveGradient(obj,w,Ix,theta,dij_interval)
            [DcIx,Dr] = getIntervalMatrices(obj,w,Ix,dij_interval);

            doseCenter = DcIx*w;
            deviation = doseCenter - obj.parameters{1};

            centerGradient = 2 * (DcIx' * deviation);
            radiusGradient = (Dr + Dr') * w;
            centerNormGradient = 2 * (DcIx' * doseCenter);

            fWGrad = obj.penalty/numel(Ix) * ...
                (centerGradient + theta * (radiusGradient - centerNormGradient));
        end

        function fFluence = bertoluzza(obj,w,Ix,theta,dij_interval)
            fFluence = computeBertoluzzaObjective(obj,w,Ix,theta,dij_interval);
        end
    end

    methods (Access = private)
        function fFluence = computeBertoluzzaObjective(obj,w,Ix,theta,dij_interval)
            [DcIx,Dr] = getIntervalMatrices(obj,w,Ix,dij_interval);

            doseCenter = DcIx*w;
            deviation = doseCenter - obj.parameters{1};

            radiusTerm = w' * Dr * w;
            centerTerm = doseCenter' * doseCenter;

            fFluence = obj.penalty/numel(Ix) * ...
                (deviation'*deviation + theta * (radiusTerm - centerTerm));
        end

        function [DcIx,Dr] = getIntervalMatrices(obj,w,Ix,dij_interval)
            matRad_cfg = MatRad_Config.instance();

            if ~isfield(dij_interval,'center') || ~isfield(dij_interval,'radius')
                matRad_cfg.dispError('dij_interval must contain center and radius fields!');
            end

            Dc = dij_interval.center;
            Dr = dij_interval.radius;
            Ix = Ix(:);

            if isempty(Ix)
                matRad_cfg.dispError('Squared Bertoluzza objective requires non-empty voxel indices!');
            end

            if size(Dc,2) ~= numel(w) || size(Dr,1) ~= numel(w) || size(Dr,2) ~= numel(w)
                matRad_cfg.dispError('dij_interval dimensions are inconsistent with the fluence vector!');
            end

            if any(Ix < 1) || any(Ix > size(Dc,1))
                matRad_cfg.dispError('Voxel indices exceed dij_interval.center dimensions!');
            end

            DcIx = Dc(Ix,:);
        end
    end

    methods (Static)
        function rob = availableRobustness()
            rob = {'INTERVAL2','INTERVAL3'};
        end
    end
end
