classdef matRad_OptimizationProblem < handle
    %matRad_OptimizationProblem Main Class for fluence optimization problems
    % Describes a standard fluence optimization problem by providing the
    % implementation of the objective & constraint function/gradient wrappers
    % and managing the mapping and backprojection of the respective dose-
    % related quantity
    %
    % References
    %   [1] https://doi.org/10.1093/imanum/draa038
    %   [2] https://doi.org/10.1002/mp.14148
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        BP
        quantityOpt = '';
        
        % Maximum approximation
        useMaxApprox = 'logsumexp'; %'pnorm'; %'cheapCOWC'; %'logsumexp'; %'none';
        
        % Parameters for pnorm approximation
        p = 30; %Can be chosen larger (closer to maximum) or smaller (closer to mean). Only tested 20 >= p >= 1
        
        % Paremeters for CheapCOWC
        p1=1;
        p2=1;
        
        % Paremeters for INTERVAL
        theta1=30;
        theta2=0.95;
        dij_interval=struct();
        cache = struct();
        
    end
    
    methods
        function obj = matRad_OptimizationProblem(backProjection)
            obj.BP = backProjection;
        end
        
        %Objective function declaration
        fVal = matRad_objectiveFunction(optiProb,w,dij,cst)

        %Objective gradient declaration
        fGrad = matRad_objectiveGradient(optiProb,w,dij,cst)
        
        %Constraint function declaration
        cVal = matRad_constraintFunctions(optiProb,w,dij,cst)
        
        %Constraint Jacobian declaration
        cJacob = matRad_constraintJacobian(optiProb,w,dij,cst)
        
        %Jacobian Structure
        jacobStruct = matRad_getJacobianStructure(optiProb,w,dij,cst)
        
        [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
        
        function lb = lowerBounds(optiProb,w)
            lb = zeros(size(w));
        end
        
        function ub = upperBounds(optiProb,w)
            ub = Inf * ones(size(w));
        end

        function [d_center, d_radius, fluenceGradient_center, fluenceGradient_radius] = getDoseInterval(optiProb, cst, structureIdx, w)
            % Check if cached values exist for the given structure index
            fieldName = sprintf('s%d', structureIdx);
            if isfield(optiProb.cache, fieldName)
                cached = optiProb.cache.(fieldName);
                d_center = cached.d_center;
                d_radius = cached.d_radius;
                fluenceGradient_center = cached.fluenceGradient_center;
                fluenceGradient_radius = cached.fluenceGradient_radius;
                return;
            end
    
            % Otherwise, compute and cache
            ixContour=1;
            subIx = cst{structureIdx,4}{ixContour};

            Dc = optiProb.dij_interval.center;
            d_center = Dc * w;
            d_center = d_center(subIx);

            fluenceGradient_center = Dc(subIx, :); 

            [~, Ix] = ismember(subIx, optiProb.dij_interval.OARSubIx);
            U = optiProb.dij_interval.U(Ix);
            S = optiProb.dij_interval.S(Ix);
            V = optiProb.dij_interval.V(Ix);
    
            nVoxels = numel(Ix);
            d_radius = zeros(nVoxels, 1);
            fluenceGradient_radius = zeros(nVoxels, numel(w));
            epsilon = 1e-8;

            if nVoxels > 500
                nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
                
                % Fallback para pruebas locales
                if isnan(nWorkers) || nWorkers == 0
                    nCores = feature('numcores');
                    nWorkers = max(1, nCores);
                end
    
                % Inicia el parpool solo si no está abierto
                if isempty(gcp('nocreate'))
                    parpool('local', nWorkers);
                end

                parfor it = 1:nVoxels
                    Dr = U{it} * S{it} * (V{it})';
                    tmp = Dr * w;
                    d_r = sqrt(w' * tmp);
                    if d_r > epsilon
                        d_radius(it) = d_r;
                        fluenceGradient_radius(it, :) = tmp' / d_r;
                    end
                end
            else
                for it = 1:nVoxels
                    Dr = U{it} * S{it} * (V{it})';
                    tmp = Dr * w;
                    d_r = sqrt(w' * tmp);
                    if d_r > epsilon
                        d_radius(it) = d_r;
                        fluenceGradient_radius(it, :) = tmp' / d_r;
                    end
                end
            end
    
            % Cache results
            optiProb.cache.w = w;
            optiProb.cache.(fieldName).d_center = d_center;
            optiProb.cache.(fieldName).d_radius = d_radius;
            optiProb.cache.(fieldName).fluenceGradient_center = fluenceGradient_center;
            optiProb.cache.(fieldName).fluenceGradient_radius = fluenceGradient_radius;
        end
    
        function clearCache(optiProb)
            % Clears the entire cache
            optiProb.cache = struct();
        end
    end
    
    methods (Access = protected)
        function [val,grad] = logSumExp(optiProb,fVals)
            % [1] stable log sum exp trick
            [fMax,ixMax] = max(fVals(:));
            
            ix = true(numel(fVals),1);
            ix(ixMax) = 0;
            
            tmp = exp(fVals - fMax);
            
            expSum = sum(tmp(ix));
            val = fMax + log1p(expSum); %log1p(x) = Matlab's numerically accurate log(1+x)
            
            grad = tmp ./ (1 + expSum);
        end
        
        function [val,grad] = pNorm(optiProb,fVals,n)
            % Implemented as proposed in [2] including a normalization for stability of the exponent.
            if nargout < 3
                n = numel(fVals);
            end
            
            p = optiProb.p;
            
            valMax = max(fVals(:));
            tmp = fVals./valMax;
            
            pNormVal = sum(tmp(:).^p)^(1/p);
            
            fac = (1/n)^(1/p);
            val = valMax*fac*pNormVal;
            
            grad = fac * (tmp ./ pNormVal).^(p-1);
        end
        
        % Cheap-minimax gradient calculation
        function [val,grad] = cheapCOWC(optiProb,fVals,fProb)
            
            val=optiProb.summaxk(fVals,fProb);
            grad=gradest(@(x) optiProb.summaxk(x,fProb),fVals);
            
        end
        
        function val = summaxk(optiProb,fVals,fProb)
            
            [~,ixKp2] = maxk(fVals(:),optiProb.p2);
            fGrad = zeros(size(fVals));
            fGrad(ixKp2) = 1;
            
            if (optiProb.p1 > 1)
                [~,ixKp1] = maxk(fVals(:),optiProb.p1-1);
                fGrad(ixKp1) = 0;
            end
            
            probSum=0;
            
            for s = 1:numel(fVals)
                if fGrad(s) ~= 0
                    probSum = probSum + fProb(s);
                end
            end
            
            val=0;
            
            for s = 1:numel(fVals)
                if fGrad(s) ~= 0
                    val = val + fProb(s)/probSum * fVals(s);
                end
            end
        end    
    end
end

