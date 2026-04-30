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
    %   [3] https://doi.org/10.1002/mp.17709
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        BP
        quantityOpt = '';
        useMaxApprox = 'logsumexp'; %'pnorm'; %'logsumexp'; %'none';
        p = 30; %Can be chosen larger (closer to maximum) or smaller (closer to mean). Only tested 20 >= p >= 1
        p1 = 1; %First Cheap-Minimax rank included in c-COWC aggregation
        p2 = 1; %Last Cheap-Minimax rank included in c-COWC aggregation
        theta1 = 30; %Bertoluzza interval target weight for INTERVAL2/INTERVAL3 objectives
        theta2 = 0.95; %OAR interval radius weight for INTERVAL3 objectives
        dij_interval = struct(); %Interval dose influence data for INTERVAL2/INTERVAL3 objectives
        intervalCache = struct(); %Cache for INTERVAL3 OAR interval dose data

        minimumW = NaN;
        maximumW = NaN;
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

        function validateCheapCOWCParameters(optiProb,numScenarios)
            matRad_cfg = MatRad_Config.instance();

            if nargin < 2 || isempty(numScenarios)
                numScenarios = numel(optiProb.BP.scenarios);
            end

            if ~isscalar(optiProb.p1) || ~isnumeric(optiProb.p1) || ...
               optiProb.p1 < 1 || fix(optiProb.p1) ~= optiProb.p1
                matRad_cfg.dispError('Cheap-Minimax parameter p1 must be a positive integer!');
            end

            if ~isscalar(optiProb.p2) || ~isnumeric(optiProb.p2) || ...
               optiProb.p2 < 1 || fix(optiProb.p2) ~= optiProb.p2
                matRad_cfg.dispError('Cheap-Minimax parameter p2 must be a positive integer!');
            end

            if optiProb.p1 > optiProb.p2
                matRad_cfg.dispError('Cheap-Minimax parameter p1 must be smaller than or equal to p2!');
            end

            if optiProb.p2 > numScenarios
                matRad_cfg.dispError('Cheap-Minimax parameter p2 exceeds the number of optimization scenarios!');
            end
        end

        function validateIntervalConfiguration(optiProb,cst,w)
            matRad_cfg = MatRad_Config.instance();
            [hasInterval,needsRadius,needsOARInterval] = optiProb.getIntervalRequirements(cst);

            if ~hasInterval
                return;
            end

            if ~isscalar(optiProb.theta1) || ~isnumeric(optiProb.theta1) || ...
               ~isfinite(optiProb.theta1) || optiProb.theta1 < 0
                matRad_cfg.dispError('INTERVAL theta1 must be a finite non-negative scalar!');
            end

            if ~isscalar(optiProb.theta2) || ~isnumeric(optiProb.theta2) || ...
               ~isfinite(optiProb.theta2) || optiProb.theta2 < 0
                matRad_cfg.dispError('INTERVAL theta2 must be a finite non-negative scalar!');
            end

            if ~isfield(optiProb.dij_interval,'center')
                matRad_cfg.dispError('INTERVAL optimization requires dij_interval.center!');
            end

            if nargin >= 3 && ~isempty(w) && size(optiProb.dij_interval.center,2) ~= numel(w)
                matRad_cfg.dispError('dij_interval.center dimensions are inconsistent with the fluence vector!');
            end

            optiProb.validateIntervalQuantity();

            if needsRadius
                optiProb.validateIntervalRadius(w);
            end

            if needsOARInterval
                optiProb.validateOARIntervalData(cst,w);
            end
        end
        
        function lb = lowerBounds(optiProb,w)
            minW = optiProb.minimumW;
            if isnan(minW)
                lb = zeros(size(w));
            elseif isscalar(minW)
                lb = minW*ones(size(w));
            elseif isvector(minW) && all(size(minW) == size(w))
                lb = minW;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Minimum Bounds for Optimization Problem could not be set!');
            end
        end
        
        function ub = upperBounds(optiProb,w)
            maxW = optiProb.maximumW;
            if isnan(maxW)
                ub = Inf(size(w));
            elseif isscalar(maxW)
                ub = maxW*ones(size(w));
            elseif isvector(maxW) && all(size(maxW) == size(w))
                ub = maxW;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Maximum Bounds for Optimization Problem could not be set!');
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

        function [val,grad] = cheapCOWC(optiProb,fVals,fProb)
            % Cheap-Minimax scenario aggregation for c-COWC [3].
            fShape = size(fVals);
            fVals = fVals(:);
            fProb = fProb(:);

            optiProb.validateCheapCOWCInputs(fVals,fProb);

            val = optiProb.summaxk(fVals,fProb);

            if nargout > 1
                matRad_cfg = MatRad_Config.instance();
                if exist('gradest','file') ~= 2
                    matRad_cfg.dispError('gradest from thirdParty/DERIVESTsuite is required for c-COWC optimization!');
                end

                grad = gradest(@(x) optiProb.summaxk(x,fProb),fVals);
                grad = reshape(grad,fShape);
            end
        end

        function val = summaxk(optiProb,fVals,fProb)
            fVals = fVals(:);
            fProb = fProb(:);

            optiProb.validateCheapCOWCInputs(fVals,fProb);

            [~,sortIx] = sort(fVals,'descend');
            selectedIx = sortIx(optiProb.p1:optiProb.p2);
            probSum = sum(fProb(selectedIx));

            if probSum <= 0
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Cheap-Minimax selected scenario ranks must have positive probability!');
            end

            val = sum(fProb(selectedIx).*fVals(selectedIx)) / probSum;
        end

        function validateCheapCOWCInputs(optiProb,fVals,fProb)
            matRad_cfg = MatRad_Config.instance();

            if numel(fVals) ~= numel(fProb)
                matRad_cfg.dispError('Cheap-Minimax values and probabilities must have equal length!');
            end

            if any(~isfinite(fVals))
                matRad_cfg.dispError('Cheap-Minimax objective values must be finite!');
            end

            if any(~isfinite(fProb)) || any(fProb < 0)
                matRad_cfg.dispError('Cheap-Minimax scenario probabilities must be finite and non-negative!');
            end

            optiProb.validateCheapCOWCParameters(numel(fVals));
        end

        function [dCenter,dRadius,fluenceGradientCenter,fluenceGradientRadius] = getOARDoseInterval(optiProb,cst,structureIdx,contourIx,w)
            optiProb.refreshIntervalCache(w);
            [subIx,refContourIx] = optiProb.getIntervalStructureVoxelIndices(cst,structureIdx);
            if ~isempty(contourIx) && contourIx ~= refContourIx
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('INTERVAL3 OAR contour index must match dij_interval.refScen!');
            end
            contourIx = refContourIx;

            fieldName = sprintf('s%d_c%d',structureIdx,contourIx);
            if isfield(optiProb.intervalCache,fieldName)
                cachedInterval = optiProb.intervalCache.(fieldName);
                dCenter = cachedInterval.dCenter;
                dRadius = cachedInterval.dRadius;
                fluenceGradientCenter = cachedInterval.fluenceGradientCenter;
                fluenceGradientRadius = cachedInterval.fluenceGradientRadius;
                return;
            end

            Dc = optiProb.dij_interval.center;
            dCenter = Dc(subIx,:)*w;
            fluenceGradientCenter = Dc(subIx,:);

            [~,intervalIx] = ismember(subIx,optiProb.dij_interval.OARSubIx(:));
            dRadius = zeros(numel(subIx),1);
            fluenceGradientRadius = zeros(numel(subIx),numel(w));
            epsilon = 1e-12;
            matRad_cfg = MatRad_Config.instance();

            for ix = 1:numel(subIx)
                U = optiProb.dij_interval.U{intervalIx(ix)};
                S = optiProb.dij_interval.S{intervalIx(ix)};
                V = optiProb.dij_interval.V{intervalIx(ix)};

                radiusGradientTmp = U*(S*(V'*w));
                radiusSquared = w'*radiusGradientTmp;

                if radiusSquared < 0 && abs(radiusSquared) < epsilon
                    radiusSquared = 0;
                elseif radiusSquared < 0
                    matRad_cfg.dispError('INTERVAL3 OAR radius matrix must be positive semidefinite!');
                end

                dRadius(ix) = sqrt(radiusSquared);

                if dRadius(ix) > epsilon
                    fluenceGradientRadius(ix,:) = radiusGradientTmp' / dRadius(ix);
                end
            end

            cachedInterval.dCenter = dCenter;
            cachedInterval.dRadius = dRadius;
            cachedInterval.fluenceGradientCenter = fluenceGradientCenter;
            cachedInterval.fluenceGradientRadius = fluenceGradientRadius;
            optiProb.intervalCache.(fieldName) = cachedInterval;
        end

        function refreshIntervalCache(optiProb,w)
            if ~isfield(optiProb.intervalCache,'w') || ~isequal(optiProb.intervalCache.w,w)
                optiProb.intervalCache = struct();
                optiProb.intervalCache.w = w;
            end
        end

        function [hasInterval,needsRadius,needsOARInterval] = getIntervalRequirements(optiProb,cst)
            matRad_cfg = MatRad_Config.instance();
            hasInterval = false;
            needsRadius = false;
            needsOARInterval = false;

            for i = 1:size(cst,1)
                if isempty(cst{i,4}) || ...
                   ~(isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET'))
                    continue;
                end

                for j = 1:numel(cst{i,6})
                    objective = cst{i,6}{j};
                    if ~isa(objective,'DoseObjectives.matRad_DoseObjective') || ...
                       ~any(strcmp(objective.robustness,{'INTERVAL2','INTERVAL3'}))
                        continue;
                    end

                    optiProb.getIntervalStructureVoxelIndices(cst,i);
                    hasInterval = true;

                    if isequal(cst{i,3},'TARGET')
                        if ~isa(objective,'DoseObjectives.matRad_SquaredBertoluzzaDeviation')
                            matRad_cfg.dispError('INTERVAL target objectives require DoseObjectives.matRad_SquaredBertoluzzaDeviation!');
                        end
                        needsRadius = true;
                    elseif strcmp(objective.robustness,'INTERVAL3')
                        needsOARInterval = true;
                    end
                end
            end
        end

        function validateIntervalRadius(optiProb,w)
            matRad_cfg = MatRad_Config.instance();

            if ~isfield(optiProb.dij_interval,'radius')
                matRad_cfg.dispError('INTERVAL target objectives require dij_interval.radius!');
            end

            if nargin >= 2 && ~isempty(w) && ...
               (size(optiProb.dij_interval.radius,1) ~= numel(w) || ...
                size(optiProb.dij_interval.radius,2) ~= numel(w))
                matRad_cfg.dispError('dij_interval.radius dimensions are inconsistent with the fluence vector!');
            end
        end

        function validateIntervalQuantity(optiProb)
            if ~isfield(optiProb.dij_interval,'quantity') || isempty(optiProb.dij_interval.quantity) || ...
               isempty(optiProb.quantityOpt)
                return;
            end

            intervalQuantity = char(optiProb.dij_interval.quantity);
            optimizationQuantity = char(optiProb.quantityOpt);

            if strcmpi(intervalQuantity,optimizationQuantity)
                return;
            end

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError(['dij_interval quantity ''%s'' is inconsistent with ' ...
                'optimization quantity ''%s''.'],intervalQuantity,optimizationQuantity);
        end

        function validateOARIntervalData(optiProb,cst,w)
            matRad_cfg = MatRad_Config.instance();
            requiredFields = {'OARSubIx','U','S','V'};

            for f = 1:numel(requiredFields)
                if ~isfield(optiProb.dij_interval,requiredFields{f})
                    matRad_cfg.dispError('INTERVAL3 OAR objectives require dij_interval.%s!',requiredFields{f});
                end
            end

            OARSubIx = optiProb.dij_interval.OARSubIx(:);
            if ~iscell(optiProb.dij_interval.U) || ~iscell(optiProb.dij_interval.S) || ...
               ~iscell(optiProb.dij_interval.V) || ...
               numel(optiProb.dij_interval.U) ~= numel(OARSubIx) || ...
               numel(optiProb.dij_interval.S) ~= numel(OARSubIx) || ...
               numel(optiProb.dij_interval.V) ~= numel(OARSubIx)
                matRad_cfg.dispError('INTERVAL3 OAR SVD fields must be cell arrays matching dij_interval.OARSubIx!');
            end

            for i = 1:size(cst,1)
                if isempty(cst{i,4}) || ~isequal(cst{i,3},'OAR')
                    continue;
                end

                hasInterval3OAR = false;
                for j = 1:numel(cst{i,6})
                    objective = cst{i,6}{j};
                    if isa(objective,'DoseObjectives.matRad_DoseObjective') && strcmp(objective.robustness,'INTERVAL3')
                        hasInterval3OAR = true;
                        break;
                    end
                end

                if ~hasInterval3OAR
                    continue;
                end

                [subIx,~] = optiProb.getIntervalStructureVoxelIndices(cst,i);
                [isCovered,intervalIx] = ismember(subIx,OARSubIx);
                if any(~isCovered)
                    matRad_cfg.dispError('INTERVAL3 OAR objectives require all OAR voxels in dij_interval.OARSubIx!');
                end

                for ix = intervalIx(:)'
                    U = optiProb.dij_interval.U{ix};
                    S = optiProb.dij_interval.S{ix};
                    V = optiProb.dij_interval.V{ix};

                    if size(U,1) ~= numel(w) || size(V,1) ~= numel(w) || ...
                       size(S,1) ~= size(S,2) || size(U,2) ~= size(S,1) || size(V,2) ~= size(S,2)
                        matRad_cfg.dispError('INTERVAL3 OAR SVD dimensions are inconsistent with the fluence vector!');
                    end
                end
            end
        end

        function [subIx,contourIx] = getIntervalStructureVoxelIndices(optiProb,cst,structureIdx)
            matRad_cfg = MatRad_Config.instance();
            contourIx = 1;

            if isfield(optiProb.dij_interval,'refScen') && ...
               ~isempty(optiProb.dij_interval.refScen)
                contourIx = optiProb.dij_interval.refScen;
            end

            if ~isnumeric(contourIx) || ~isscalar(contourIx) || ...
               ~isfinite(contourIx) || contourIx < 1 || fix(contourIx) ~= contourIx
                matRad_cfg.dispError('dij_interval.refScen must be a positive integer scalar!');
            end
            contourIx = double(contourIx);

            if size(cst,2) < 4 || structureIdx > size(cst,1) || ...
               isempty(cst{structureIdx,4}) || ~iscell(cst{structureIdx,4}) || ...
               numel(cst{structureIdx,4}) < contourIx || ...
               isempty(cst{structureIdx,4}{contourIx})
                matRad_cfg.dispError(['INTERVAL objectives require non-empty structure ', ...
                    'voxel indices for reference CT scenario %d.'],contourIx);
            end

            subIx = cst{structureIdx,4}{contourIx}(:);
        end
    end                
end
