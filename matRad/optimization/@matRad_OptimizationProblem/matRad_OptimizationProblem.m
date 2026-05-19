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
    % Copyright 2020-2026 the matRad development team.
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
        p1 = 1; % First Cheap-Minimax rank included in c-COWC aggregation
        p2 = 1; % Last Cheap-Minimax rank included in c-COWC aggregation
        theta1 = 30; %Bertoluzza interval target weight for INTERVAL2/INTERVAL3 objectives
        theta2 = 1.0; %OAR interval radius weight for INTERVAL3 objectives
        dij_interval = struct(); %Interval dose influence data for INTERVAL2/INTERVAL3 objectives
        intervalCache = struct(); %Cache for INTERVAL3 OAR interval dose data
        dij_prob = struct(); %Scenario-free probabilistic data for PROB/PROB2 optimization
        prob2Cache = struct(); %Cache for PROB/PROB2 expected dose and variance data

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

        function validateProb2Configuration(optiProb,cst,w)
            [hasProb2,needsOmega] = optiProb.getProb2Requirements(cst);

            if ~hasProb2
                return;
            end

            optiProb.validateProb2Data(w);

            if needsOmega
                for i = 1:size(cst,1)
                    if optiProb.structureNeedsProb2Omega(cst,i)
                        optiProb.validateProb2OmegaData(cst,i,w);
                    end
                end
            end
        end

        function stats = GetResultProbabilistic(optiProb,w,dij,cst,structureIdx) %#ok<INUSD>
            stats = optiProb.getProb2DoseStats(w,cst,structureIdx);
        end

        function stats = GetResultInterval(optiProb,w,cst,structureIdx,objective)
            stats = optiProb.getIntervalDoseStats(w,cst,structureIdx,objective);
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

    methods (Access = private)
        function refreshIntervalCache(optiProb,w)
            if ~isfield(optiProb.intervalCache,'w') || ~isequal(optiProb.intervalCache.w,w)
                optiProb.intervalCache = struct();
                optiProb.intervalCache.w = w;
            end
        end

        function refreshProb2Cache(optiProb,w)
            if ~isfield(optiProb.prob2Cache,'w') || ~isequal(optiProb.prob2Cache.w,w)
                optiProb.prob2Cache = struct();
                optiProb.prob2Cache.w = w;
            end
        end

        function stats = getProb2DoseStats(optiProb,w,cst,structureIdx)
            optiProb.refreshProb2Cache(w);
            [subIx,contourIx] = optiProb.getProb2StructureVoxelIndices(cst,structureIdx);
            fieldName = sprintf('s%d_c%d',structureIdx,contourIx);

            if isfield(optiProb.prob2Cache,fieldName)
                stats = optiProb.prob2Cache.(fieldName);
                return;
            end

            optiProb.validateProb2Data(w);
            expectedRows = optiProb.dij_prob.expected(subIx,:);

            stats = struct();
            stats.subIx = subIx;
            stats.contourIx = contourIx;
            stats.expectedRows = expectedRows;
            stats.dExp = expectedRows * w;
            stats.Omega = optiProb.getProb2OmegaForStructure(cst,structureIdx,subIx);
            stats.omegaW = [];
            stats.meanVariance = [];
            stats.gradMeanVariance = [];

            if ~isempty(stats.Omega)
                stats.omegaW = stats.Omega * w;
                stats.meanVariance = full(w' * stats.omegaW) / numel(subIx);
                stats.meanVariance = optiProb.validateProb2MeanVariance( ...
                    stats.meanVariance,w,stats.omegaW);
                stats.gradMeanVariance = 2 .* stats.omegaW ./ numel(subIx);
            end

            optiProb.prob2Cache.(fieldName) = stats;
        end

        function stats = getIntervalDoseStats(optiProb,w,cst,structureIdx,objective,contourIx)
            if nargin < 6
                contourIx = [];
            end

            matRad_cfg = MatRad_Config.instance();
            robustness = objective.robustness;
            [subIx,refContourIx] = optiProb.getIntervalStructureVoxelIndices(cst,structureIdx);
            if ~isempty(contourIx) && contourIx ~= refContourIx
                matRad_cfg.dispError('INTERVAL contour index must match dij_interval.refScen!');
            end

            optiProb.validateIntervalQuantity();
            if isequal(cst{structureIdx,3},'TARGET')
                if ~isa(objective,'DoseObjectives.matRad_SquaredBertoluzzaDeviation')
                    matRad_cfg.dispError('INTERVAL target objectives require DoseObjectives.matRad_SquaredBertoluzzaDeviation!');
                end
                optiProb.validateIntervalRadius(w);
            elseif strcmp(robustness,'INTERVAL3')
                optiProb.validateOARIntervalData(cst,w);
            end

            stats = struct();
            stats.subIx = subIx;
            stats.contourIx = refContourIx;
            stats.centerRows = optiProb.dij_interval.center(subIx,:);
            stats.centerDose = stats.centerRows * w;
            stats.radiusDose = [];
            stats.gradRadius = [];
            stats.doseForObjective = stats.centerDose;
            stats.gradDoseForObjective = stats.centerRows;
            stats.radiusMatrix = [];

            if isequal(cst{structureIdx,3},'TARGET')
                if isfield(optiProb.dij_interval,'radius')
                    stats.radiusMatrix = optiProb.dij_interval.radius;
                end
                return;
            end

            if strcmp(robustness,'INTERVAL2')
                return;
            end

            optiProb.refreshIntervalCache(w);
            fieldName = sprintf('s%d_c%d_interval3',structureIdx,refContourIx);
            if isfield(optiProb.intervalCache,fieldName)
                cachedInterval = optiProb.intervalCache.(fieldName);
                stats.radiusDose = cachedInterval.radiusDose;
                stats.gradRadius = cachedInterval.gradRadius;
                stats.doseForObjective = stats.centerDose + ...
                    optiProb.theta2 .* stats.radiusDose;
                stats.gradDoseForObjective = stats.centerRows + ...
                    optiProb.theta2 .* stats.gradRadius;
                return;
            end

            [~,intervalIx] = ismember(subIx,optiProb.dij_interval.OARSubIx(:));
            if any(intervalIx == 0)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('INTERVAL3 OAR objectives require all OAR voxels in dij_interval.OARSubIx!');
            end

            radiusDose = zeros(numel(subIx),1);
            gradRadius = zeros(numel(subIx),numel(w));
            epsilon = 1e-12;

            for ix = 1:numel(subIx)
                radiusFactor = optiProb.dij_interval.OARRadiusFactor{intervalIx(ix)};
                projectedRadius = radiusFactor' * w;
                radiusDose(ix) = norm(projectedRadius);

                if radiusDose(ix) > epsilon
                    gradRadius(ix,:) = (radiusFactor * projectedRadius)' ./ radiusDose(ix);
                end
            end

            stats.radiusDose = radiusDose;
            stats.gradRadius = gradRadius;
            stats.doseForObjective = stats.centerDose + optiProb.theta2 .* stats.radiusDose;
            stats.gradDoseForObjective = stats.centerRows + optiProb.theta2 .* stats.gradRadius;

            cachedInterval.radiusDose = stats.radiusDose;
            cachedInterval.gradRadius = stats.gradRadius;
            optiProb.intervalCache.(fieldName) = cachedInterval;
        end

        function validateProb2Data(optiProb,w)
            matRad_cfg = MatRad_Config.instance();

            if ~isfield(optiProb.dij_prob,'expected')
                matRad_cfg.dispError('PROB/PROB2 optimization requires dij_prob.expected!');
            end

            if nargin >= 2 && ~isempty(w) && ...
                    size(optiProb.dij_prob.expected,2) ~= numel(w)
                matRad_cfg.dispError('dij_prob.expected dimensions are inconsistent with the fluence vector!');
            end

            optiProb.validateProb2Quantity();
        end

        function validateProb2Quantity(optiProb)
            if ~optiProb.hasInconsistentQuantity(optiProb.dij_prob)
                return;
            end

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError(['dij_prob quantity ''%s'' is inconsistent with ', ...
                                  'optimization quantity ''%s''.'], ...
                                  optiProb.dij_prob.quantity,optiProb.quantityOpt);
        end

        function validateProb2OmegaData(optiProb,cst,structureIdx,w)
            [subIx,~] = optiProb.getProb2StructureVoxelIndices(cst,structureIdx);
            Omega = optiProb.getProb2OmegaForStructure(cst,structureIdx,subIx);
            matRad_cfg = MatRad_Config.instance();

            if isempty(Omega)
                matRad_cfg.dispError('Probabilistic mean variance requires dij_prob.Omega{%d}!',structureIdx);
            end

            if nargin >= 4 && ~isempty(w) && ...
                    (size(Omega,1) ~= numel(w) || size(Omega,2) ~= numel(w))
                matRad_cfg.dispError('dij_prob.Omega{%d} dimensions are inconsistent with the fluence vector!',structureIdx);
            end
        end

        function meanVariance = validateProb2MeanVariance(~,meanVariance,w,omegaW)
            epsilon = 1e-10 * max(1,norm(w) * norm(omegaW));
            if meanVariance < 0 && abs(meanVariance) <= epsilon
                meanVariance = 0;
            elseif meanVariance < 0
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('dij_prob.Omega matrix must be positive semidefinite for the current fluence vector!');
            end
        end

        function Omega = getProb2OmegaForStructure(optiProb,cst,structureIdx,subIx)
            Omega = [];
            if ~isfield(optiProb.dij_prob,'Omega') || ...
                    numel(optiProb.dij_prob.Omega) < structureIdx || ...
                    isempty(optiProb.dij_prob.Omega{structureIdx})
                return;
            end

            if isfield(optiProb.dij_prob,'voiSubIx') && ...
                    numel(optiProb.dij_prob.voiSubIx) >= structureIdx && ...
                    ~isempty(optiProb.dij_prob.voiSubIx{structureIdx}) && ...
                    ~isequal(optiProb.dij_prob.voiSubIx{structureIdx}(:),subIx(:))
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('dij_prob.voiSubIx{%d} does not match the selected cst voxels.',structureIdx);
            end

            if isempty(cst{structureIdx,4})
                return;
            end
            Omega = optiProb.dij_prob.Omega{structureIdx};
        end

        function [hasProb2,needsOmega] = getProb2Requirements(optiProb,cst)
            matRad_cfg = MatRad_Config.instance();
            hasProb2 = false;
            needsOmega = false;

            for i = 1:size(cst,1)
                if isempty(cst{i,4}) || ...
                        ~(isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET'))
                    continue;
                end

                for j = 1:numel(cst{i,6})
                    optiFunc = cst{i,6}{j};
                    if ~isa(optiFunc,'matRad_DoseOptimizationFunction') || ...
                            ~any(strcmp(optiFunc.robustness,{'PROB','PROB2'}))
                        continue;
                    end

                    if strcmp(optiFunc.robustness,'PROB') && ...
                            isa(optiFunc,'DoseConstraints.matRad_MinMaxMeanVariance')
                        matRad_cfg.dispError('MinMaxMeanVariance constraints are only supported for PROB2 robustness!');
                    end

                    optiProb.getProb2StructureVoxelIndices(cst,i);
                    hasProb2 = true;

                    if strcmp(optiFunc.robustness,'PROB') && ...
                            isa(optiFunc,'DoseObjectives.matRad_DoseObjective')
                        needsOmega = true;
                    elseif strcmp(optiFunc.robustness,'PROB2') && ...
                            (isa(optiFunc,'DoseObjectives.matRad_MeanVariance') || ...
                            isa(optiFunc,'DoseConstraints.matRad_MinMaxMeanVariance'))
                        needsOmega = true;
                    end
                end
            end
        end

        function needsOmega = structureNeedsProb2Omega(~,cst,structureIdx)
            needsOmega = false;

            if isempty(cst{structureIdx,4}) || ...
                    ~(isequal(cst{structureIdx,3},'OAR') || ...
                    isequal(cst{structureIdx,3},'TARGET'))
                return;
            end

            for j = 1:numel(cst{structureIdx,6})
                optiFunc = cst{structureIdx,6}{j};
                if ~isa(optiFunc,'matRad_DoseOptimizationFunction')
                    continue;
                end

                if strcmp(optiFunc.robustness,'PROB') && ...
                        isa(optiFunc,'DoseObjectives.matRad_DoseObjective')
                    needsOmega = true;
                    return;
                end

                if strcmp(optiFunc.robustness,'PROB2') && ...
                        (isa(optiFunc,'DoseObjectives.matRad_MeanVariance') || ...
                        isa(optiFunc,'DoseConstraints.matRad_MinMaxMeanVariance'))
                    needsOmega = true;
                    return;
                end
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
            if ~optiProb.hasInconsistentQuantity(optiProb.dij_interval)
                return;
            end

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError(['dij_interval quantity ''%s'' is inconsistent with ', ...
                                  'optimization quantity ''%s''.'], ...
                                  optiProb.dij_interval.quantity,optiProb.quantityOpt);
        end

        function validateOARIntervalData(optiProb,cst,w)
            matRad_cfg = MatRad_Config.instance();
            requiredFields = {'OARSubIx','OARRadiusFactor','OARRadiusRank'};

            for f = 1:numel(requiredFields)
                if ~isfield(optiProb.dij_interval,requiredFields{f})
                    matRad_cfg.dispError('INTERVAL3 OAR objectives require dij_interval.%s!',requiredFields{f});
                end
            end

            OARSubIx = optiProb.dij_interval.OARSubIx(:);
            if ~iscell(optiProb.dij_interval.OARRadiusFactor) || ...
                    numel(optiProb.dij_interval.OARRadiusFactor) ~= numel(OARSubIx) || ...
                    numel(optiProb.dij_interval.OARRadiusRank) ~= numel(OARSubIx)
                matRad_cfg.dispError(['INTERVAL3 OAR radius fields must match ', ...
                                      'dij_interval.OARSubIx!']);
            end

            for i = 1:size(cst,1)
                if isempty(cst{i,4}) || ~isequal(cst{i,3},'OAR')
                    continue;
                end

                hasInterval3OAR = false;
                for j = 1:numel(cst{i,6})
                    objective = cst{i,6}{j};
                    if isa(objective,'DoseObjectives.matRad_DoseObjective') && ...
                            strcmp(objective.robustness,'INTERVAL3')
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
                    radiusFactor = optiProb.dij_interval.OARRadiusFactor{ix};
                    radiusRank = optiProb.dij_interval.OARRadiusRank(ix);

                    if ~isnumeric(radiusRank) || ~isscalar(radiusRank) || ...
                            ~isfinite(radiusRank) || radiusRank < 0 || ...
                            fix(radiusRank) ~= radiusRank
                        matRad_cfg.dispError('INTERVAL3 OAR radius rank must be a non-negative integer!');
                    end

                    if ~isnumeric(radiusFactor) || size(radiusFactor,1) ~= numel(w) || ...
                            size(radiusFactor,2) ~= radiusRank
                        matRad_cfg.dispError(['INTERVAL3 OAR radius factor dimensions ', ...
                                              'are inconsistent with the fluence vector!']);
                    end

                    if any(~isfinite(nonzeros(radiusFactor)))
                        matRad_cfg.dispError('INTERVAL3 OAR radius factor contains non-finite values!');
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

            contourIx = optiProb.validateReferenceScenario(contourIx,'dij_interval.refScen');

            if size(cst,2) < 4 || structureIdx > size(cst,1) || ...
                    isempty(cst{structureIdx,4}) || ~iscell(cst{structureIdx,4}) || ...
                    numel(cst{structureIdx,4}) < contourIx || ...
                    isempty(cst{structureIdx,4}{contourIx})
                matRad_cfg.dispError(['INTERVAL objectives require non-empty structure ', ...
                                      'voxel indices for reference CT scenario %d.'],contourIx);
            end

            subIx = cst{structureIdx,4}{contourIx}(:);
        end

        function [subIx,contourIx] = getProb2StructureVoxelIndices(optiProb,cst,structureIdx)
            matRad_cfg = MatRad_Config.instance();
            contourIx = 1;

            if isfield(optiProb.dij_prob,'refScen') && ...
                    ~isempty(optiProb.dij_prob.refScen)
                contourIx = optiProb.dij_prob.refScen;
            end

            contourIx = optiProb.validateReferenceScenario(contourIx,'dij_prob.refScen');

            if size(cst,2) < 4 || structureIdx > size(cst,1) || ...
                    isempty(cst{structureIdx,4}) || ~iscell(cst{structureIdx,4}) || ...
                    numel(cst{structureIdx,4}) < contourIx || ...
                    isempty(cst{structureIdx,4}{contourIx})
                matRad_cfg.dispError(['Probabilistic optimization requires non-empty structure ', ...
                                      'voxel indices for reference CT scenario %d.'],contourIx);
            end

            subIx = cst{structureIdx,4}{contourIx}(:);
        end

        function contourIx = validateReferenceScenario(~,contourIx,fieldName)
            if ~isnumeric(contourIx) || ~isscalar(contourIx) || ...
                    ~isfinite(contourIx) || contourIx < 1 || fix(contourIx) ~= contourIx
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('%s must be a positive integer scalar!',fieldName);
            end
            contourIx = double(contourIx);
        end

        function tf = hasInconsistentQuantity(optiProb,payload)
            tf = false;
            if ~isfield(payload,'quantity') || isempty(payload.quantity) || isempty(optiProb.quantityOpt)
                return;
            end

            payloadQuantity = payload.quantity;
            if isfield(payload,'optimizationQuantity') && ~isempty(payload.optimizationQuantity)
                payloadQuantity = payload.optimizationQuantity;
            end

            tf = ~strcmpi(optiProb.canonicalQuantity(payloadQuantity), ...
                          optiProb.canonicalQuantity(optiProb.quantityOpt));
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

            val = sum(fProb(selectedIx) .* fVals(selectedIx)) / probSum;
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
    end                

    methods (Static, Access = private)
        function quantity = canonicalQuantity(quantity)
            if isstring(quantity) && isscalar(quantity)
                quantity = char(quantity);
            end
        end
    end
end
