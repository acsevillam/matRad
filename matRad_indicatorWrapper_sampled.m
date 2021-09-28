function [dvh,qi] = matRad_indicatorWrapper_sampled(cst,pln,resultGUI,refGy,refVol,param)
% matRad indictor wrapper
% 
% call
%   matRad_calcIndicators(cst,pln,cube,dvhType,param,refGy,refVol,lineStyleIndicator)
%
% input
%   cst:                  matRad cst struct
%   pln:                  matRad pln struct
%   resultGUI:            matRad resultGUI struct
%   refGy: (optional)     array of dose values used for V_XGy calculation
%                         default is [40 50 60]
%   refVol:(optional)     array of volumes (0-100) used for D_X calculation
%                         default is [2 5 95 98]
%                         NOTE: Call either both or none!
%
% output
%   dvh: matRad dvh result struct
%   qi:  matRad quality indicator result struct
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

if ~exist('refVol', 'var') 
    refVol = [];
end

if ~exist('refGy', 'var')
	refGy = [];
end

if exist('param','var')
	if ~isfield(param,'logLevel')
        param.logLevel = 1;
	end
else
	param.logLevel = 1;
end

dvh = cell(1,length(3));
qi = cell(1,length(3));

doseCube = resultGUI.([pln.bioParam.quantityVis]);
dvh{1,1} = matRad_calcDVH(cst,doseCube,'cum');
qi{1,1}  = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol);

doseCube = resultGUI.([pln.bioParam.quantityVis '_lower']);
dvh{1,2} = matRad_calcDVH(cst,doseCube,'cum');
qi{1,2}  = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol);

doseCube = resultGUI.([pln.bioParam.quantityVis '_upper']);
dvh{1,3} = matRad_calcDVH(cst,doseCube,'cum');
qi{1,3}  = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol);

x0=10;
y0=10;
width=600;
height=400;
figure,set(gcf,'Color',[1 1 1],'position',[x0,y0,width,height]);
%subplot(2,1,1)
matRad_showDVH_sampled(dvh,cst,pln);
%subplot(2,1,2)
%matRad_showQualityIndicators(qi{1,1});


