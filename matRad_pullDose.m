function [cst,flag] = matRad_pullDose(cst,pullingStep)
% matRad add margin function
%
% call
%   mVOIEnlarged = matRad_addMargin(mVOI,cst,vResolution,vMargin,bDiaElem)
%
% input
%   mVOI:           image stack in dimensions of X x Y x Z holding ones for
%                   object and zeros otherwise
%   cst:            matRad cst struct
%   vResolution     ct resolution
%   vMargin:        margin in mm
%   bDiaElem        if true 26-connectivity is used otherwise 6-connectivity
%
% output
%   mVOIEnlarged:   enlarged VOI
%
% References
%   -
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

flag=false;

% compute for every VOI.
for  i = 1:size(cst,1)
    if ~isempty(cst{i,4}{1}) % && isequal(cst{i,3},'OAR')
        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            % Get current objective
            objective = cst{i,6}{j};
            if(objective.dosePulling)
                if(objective.pullingStep==pullingStep)
                    % Dose pulling
                    for k = 1:numel(objective.parameters)
                        if ~isempty(objective.objectivePullingRate{k}) && ((objective.objectivePullingRate{k}>0 && objective.parameters{k}>=0) || (objective.objectivePullingRate{k}<0 && objective.parameters{k}>0) )
                            objective.parameters{k}=objective.parameters{k}+objective.objectivePullingRate{k};
                            flag=true;
                        end
                    end
        
                    if ~isempty(objective.penaltyPullingRate) && objective.penaltyPullingRate~=0
                        objective.penalty=objective.penalty+objective.penaltyPullingRate;
                        flag=true;
                    end
                    
                    cst{i,6}{j}=objective;
                end
            end

        end
    end

end


end

