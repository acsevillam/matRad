function [dvh] = matRad_multScenDVHWrapper(cst,pln,resultGUI_multScen,param)
% matRad multi scenario DVH wrapper
% 
% call
%   matRad_MultScenDVHWrapper(cst,pln,resultGUI,param)
%
% input
%   cst:                  matRad cst struct
%   pln:                  matRad pln struct
%   resultGUI_multScen:   matRad resultGUI structures for multiple
%   scenarios
%
% output
%   dvh: matRad dvh result struct
%   graphical display of all results
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('param','var')
   if ~isfield(param,'logLevel')
      param.logLevel = 1;
   end
else
   param.logLevel = 1;
end

numScen=size(resultGUI_multScen);

dvh = cell(numScen(1),numScen(2),numScen(3));

for ctScen=1:numScen(1)
    for shiftScen=1:numScen(2)
        for shiftRangeScen=1:numScen(3)  
            if(isempty(resultGUI_multScen{ctScen,shiftScen,shiftRangeScen})==false)
                doseCube = resultGUI_multScen{ctScen,shiftScen,shiftRangeScen}.([pln.bioParam.quantityOpt]);
                metadata.ctScen=ctScen;
                dvh{ctScen,shiftScen,shiftRangeScen} = matRad_calcDVH(cst,doseCube,'cum',[],metadata);
            end
        end
    end
end

figure,set(gcf,'Color',[1 1 1]);
matRad_showMultScenDVH(dvh,cst,pln);



