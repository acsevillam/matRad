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
            % Retrieves or computes the expected dose (d_center), dose uncertainty (d_radius),
            % and their respective fluence gradients for a given structure.
            %
            % Caching behavior:
            % - Returns cached d_center and fluenceGradient_center if available.
            % - Returns cached d_radius and fluenceGradient_radius if available and valid.
            % - Otherwise, computes missing quantities and stores them in the cache.
            % - Stores the current w as w_last_radius only when d_radius is computed.
            % - This function does not evaluate whether the cache is outdated;
            %   call optiProb.clearCache(w) beforehand to invalidate outdated entries.

            fieldName = sprintf('s%d', structureIdx);
        
            % Initialize cache field if not present
            if ~isfield(optiProb.cache, fieldName)
                optiProb.cache.(fieldName) = struct();
            end
        
            cacheStruct = optiProb.cache.(fieldName);
        
            % Load subIx
            ixContour = 1;
            subIx = cst{structureIdx,4}{ixContour};
        
            % Center dose (d_center) and its gradient
            if isfield(cacheStruct, 'd_center') && isfield(cacheStruct, 'fluenceGradient_center')
                d_center = cacheStruct.d_center;
                fluenceGradient_center = cacheStruct.fluenceGradient_center;
            else
                Dc = optiProb.dij_interval.center;
                d_center = Dc * w;
                d_center = d_center(subIx);
                fluenceGradient_center = Dc(subIx, :);
        
                optiProb.cache.(fieldName).d_center = d_center;
                optiProb.cache.(fieldName).fluenceGradient_center = fluenceGradient_center;
            end
        
            % Radius dose (d_radius) and its gradient
            if isfield(cacheStruct, 'd_radius') && isfield(cacheStruct, 'fluenceGradient_radius')
                d_radius = cacheStruct.d_radius;
                fluenceGradient_radius = cacheStruct.fluenceGradient_radius;
            else
                % Load SVD components
                [~, Ix] = ismember(subIx, optiProb.dij_interval.OARSubIx);
                U = optiProb.dij_interval.U(Ix);
                S = optiProb.dij_interval.S(Ix);
                V = optiProb.dij_interval.V(Ix);
        
                % Allocate memory
                nVoxels = numel(Ix);
                d_radius = zeros(nVoxels, 1);
                fluenceGradient_radius = zeros(nVoxels, numel(w));
                epsilon = 1e-8;
        
                % Parallel or serial computation
                if nVoxels > 500
                    nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
                    if isnan(nWorkers) || nWorkers == 0
                        nCores = feature('numcores');
                        nWorkers = max(1, nCores)-2;
                    end
                    if isempty(gcp('nocreate'))
                        parpool('local', nWorkers);
                    end
        
                    parfor it = 1:nVoxels
                        Dr = U{it} * S{it} * V{it}';
                        tmp = Dr * w;
                        d_r = sqrt(w' * tmp);
                        if d_r > epsilon
                            d_radius(it) = d_r;
                            fluenceGradient_radius(it, :) = tmp' / d_r;
                        end
                    end
                else
                    for it = 1:nVoxels
                        Dr = U{it} * S{it} * V{it}';
                        tmp = Dr * w;
                        d_r = sqrt(w' * tmp);
                        if d_r > epsilon
                            d_radius(it) = d_r;
                            fluenceGradient_radius(it, :) = tmp' / d_r;
                        end
                    end
                end
        
                % Store results and update last radius w
                optiProb.cache.(fieldName).d_radius = d_radius;
                optiProb.cache.(fieldName).fluenceGradient_radius = fluenceGradient_radius;
                optiProb.cache.(fieldName).w_last_radius = w;
            end
        end

        function clearCache(optiProb, w)
            % Clears cached data depending on how much 'w' has changed.
            % - Clears d_center and fluenceGradient_center if w has changed.
            % - Clears d_radius aand fluenceGradient_radius only if w differs significantly from the previous w_last_radius that generated it.
            % - Completely clears cache every 5 changes of w.
        
            % If no cache yet, create new cache
            if ~isfield(optiProb.cache, 'w')
                fullClear();
                return;
            end
        
            w_cached = optiProb.cache.w;
        
            % === CASE 1: w is identical → do nothing ===
            if isequal(w_cached, w)
                return;
            end
        
            w_last_radius_cached = optiProb.cache.w_last_radius;
            epsilon = 0.1;  % Relative difference threshold for partial clear
            tolerance = 0.0;

            % Compute relative change
            rel_diff = abs(w - w_last_radius_cached) ./ max(abs(w_last_radius_cached), eps);
        
            if sum(rel_diff > epsilon) / numel(rel_diff) <= tolerance
                % === CASE 2: (relative change <= epsilon) → clear only d_center and fluenceGradient_center ===
        
                if mod(optiProb.cache.clearCalls, 10) == 0
                    fullClear();
                else
                    clearCenterOnly();
                    optiProb.cache.clearCalls = optiProb.cache.clearCalls + 1;
                end

            else
                % === CASE 3: (relative change > epsilon) → full clear ===    
                fullClear();
            end
        
            % ===============================
            function fullClear()
                optiProb.cache = struct();
                optiProb.cache.w = w;
                optiProb.cache.w_last_radius = w;
                optiProb.cache.clearCalls = 1;
                fprintf('[interval] Full cache cleared.\n');
            end
        
            function clearCenterOnly()
                fieldList = fieldnames(optiProb.cache);
                for i = 1:numel(fieldList)
                    key = fieldList{i};
                    if ismember(key, {'w', 'clearCalls'})
                        continue;
                    end
                    structureCache = optiProb.cache.(key);
                    structureCache = rmfield_safe(structureCache, 'd_center');
                    structureCache = rmfield_safe(structureCache, 'fluenceGradient_center');
                    optiProb.cache.(key) = structureCache;
                end
            end
        
            function s = rmfield_safe(s, f)
                if isfield(s, f)
                    s = rmfield(s, f);
                end
            end
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

