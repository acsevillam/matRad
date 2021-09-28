classdef matRad_SquaredBertoluzzaDeviation < DoseObjectives.matRad_DoseObjective
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
        parameterNames = {'d^{ref}','theta'};
        parameterTypes = {'dose','theta'};
    end
    
    properties
        parameters = {60,0.1};
        penalty = 1;
    end

    
    methods
        function obj = matRad_SquaredBertoluzzaDeviation(penalty,dRef,theta)
            
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
                if nargin == 3 && isscalar(theta)
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
        function fDose = computeDoseObjectiveFunction(obj,dose)
            dose_center=dose.dose_center;
            dose_radius=dose.dose_radius;
            % deviation : dose minus prefered dose
            deviation = dose_center - obj.parameters{1};
            % claculate objective function
            fDose = obj.penalty/numel(dose_center) * ((deviation'*deviation) + ...
                obj.parameters{2} * (dose_radius'*dose_radius));
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            dose_center=dose.dose_center;
            dose_radius=dose.dose_radius;
            % deviation : Dose minus prefered dose
            deviation = dose_center - obj.parameters{1};
            % calculate delta
            fDoseGrad = obj.penalty/numel(dose_center) * (2 * deviation + ...
                2 * obj.parameters{2} * dose_radius);
        end
    end
    
    methods (Static)
        function rob = availableRobustness()
            rob = DoseObjectives.matRad_DoseObjective.availableRobustness();
            rob{end+1} = 'INTERVAL1'; %By default, no robustness is available
            rob{end+1} = 'INTERVAL2'; %By default, no robustness is available
            rob{end+1} = 'INTERVAL3'; %By default, no robustness is available
            rob{end+1} = 'INTERVAL4'; %By default, no robustness is available
        end
    end
    
end
