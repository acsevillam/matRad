classdef matRad_SquaredBertoluzzaDeviation2 < DoseObjectives.matRad_DoseObjective
    % matRad_SquaredDeviation Implements a penalized least squares objective
    %   See matRad_DoseObjective for interface description
    %
    % References
    %     -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2015 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
        name = 'Squared Deviation';
        parameterNames = {'d^{ref}','theta','dij_interval'};
        parameterTypes = {'dose','theta','dij_interval'};
    end
    
    properties
        parameters = {60,0.1};
        penalty = 1;
    end
    
    methods
        function obj = matRad_SquaredBertoluzzaDeviation2(penalty,dRef,theta,dij_interval)
            
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
                if nargin == 4 && isstruct(dij_interval)
                    obj.parameters{3} = dij_interval;
                end
                if nargin >= 3 && isscalar(theta)
                    obj.parameters{2} = theta;
                end
                
                if nargin >= 2 && isscalar(dRef)
                    obj.parameters{1} = dRef;
                end
                
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
            end
        end
        
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,w,Ix)
            Dc = obj.parameters{3}.center;
            Dr = obj.parameters{3}.radius;
            dose_center=Dc*w;
            dose_center=dose_center(Ix);

            % deviation : dose minus prefered dose
            deviation = dose_center - obj.parameters{1};
            
            % radius dose first term
            dose_radius_1 = w'*Dr*w;

            % calculate objective function
            fDose = obj.penalty/numel(dose_center) * (deviation'*deviation + obj.parameters{2} * (dose_radius_1 - dose_center'*dose_center));
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,w,Ix)
            Dc = obj.parameters{3}.center;
            Dr = obj.parameters{3}.radius;
            dose_center=Dc*w;
            dose_center=dose_center(Ix);

            % deviation : dose minus prefered dose
            deviation = dose_center - obj.parameters{1};

            % radius dose gradient first term
            dose_radius_grad_1 = Dc*(w'*Dr)';
            dose_radius_grad_1=dose_radius_grad_1(Ix);

            % calculate delta
            fDoseGrad = obj.penalty/numel(dose_center) * (2 * deviation + ...
                2 * obj.parameters{2} * (dose_radius_grad_1 - dose_center));

        end

    end
    
    methods (Static)
        function rob = availableRobustness()
            rob = DoseObjectives.matRad_DoseObjective.availableRobustness();
            rob{end+1} = 'INTERVAL2'; %By default, no robustness is available
        end
    end
    
end