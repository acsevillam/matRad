function jacobStruct = matRad_getJacobianStructure(optiProb,w,dij,cst)	
% matRad IPOPT callback: jacobian structure function for inverse planning 
% supporting max dose constraint, min dose constraint, min mean dose constraint, 
% max mean dose constraint, min EUD constraint, max EUD constraint, max DVH 
% constraint, min DVH constraint 	
% 	
% call:	
%   jacobStruct = matRad_getJacobStruct(optiProb,w,dij,cst)	
%	
% input:	
%   optiProb: matRad optimization problem
%   w:        beamlet/ pencil beam weight vector
%   dij: dose influence matrix	
%   cst: matRad cst struct	
%	
% output:	
%   jacobStruct: jacobian of constraint function	
%	
% References	
%	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016-2026 the matRad development team.
% 	
% This file is part of the matRad project. It is subject to the license 	
% terms in the LICENSE file found in the top-level directory of this 	
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 	
% of the matRad project, including this file, may be copied, modified, 	
% propagated, or distributed except according to the terms contained in the 	
% LICENSE file.	
%	
jacobStruct = sparse(0,dij.totalNumOfBixels);
matRad_cfg = MatRad_Config.instance();

tmp = false(size(dij.physicalDose{1},1),1);
for i = 1:size(cst,1)
    if ~isempty(cst{i,4}) && any(strcmp(cst{i,3},{'OAR','TARGET','EXTERNAL'}))
        for j = 1:numel(cst{i,6})
            obj = cst{i,6}{j};
            if ~isa(obj,'DoseConstraints.matRad_DoseConstraint')
                continue;
            end

            if any(strcmp(obj.robustness, {'PROB','PROB2'}))
                stats = optiProb.GetResultProbabilistic(w,dij,cst,i);
                if isa(obj,'DoseConstraints.matRad_MinMaxMeanVariance')
                    if strcmp(obj.robustness,'PROB')
                        matRad_cfg.dispError('MinMaxMeanVariance constraints are only supported for PROB2 robustness!');
                    end
                    if isempty(stats.Omega)
                        matRad_cfg.dispError('PROB2 mean variance requires dij_prob.Omega{%d}!',i);
                    end
                    jacobStruct = [jacobStruct; ...
                                   spones(sum(spones(stats.Omega),1))];
                else
                    jacobDoseStruct = obj.getDoseConstraintJacobianStructure(numel(stats.subIx));
                    jacobStruct = [jacobStruct; ...
                                   spones(jacobDoseStruct' * stats.expectedRows)];
                end
                continue;
            end

            if ~iscell(cst{i,4}) || isempty(cst{i,4}{1})
                continue;
            end

            tmp(:) = false;
            tmp(cst{i,4}{1}) = true;

            jacobDoseStruct = obj.getDoseConstraintJacobianStructure(numel(cst{i,4}{1}));
            nRows = size(jacobDoseStruct,2);
            jacobStruct = [jacobStruct; ...
                           repmat(spones(double(tmp') * dij.physicalDose{1}),nRows,1)];
        end
    end
end
