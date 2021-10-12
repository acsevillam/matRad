function matRad_showUVH(dvh,uvh,cst,pln,doseWindow,lineStyleIndicator)
% matRad uvh visualizaion
% 
% call
%   matRad_showUVH(uvh,cst)
%   matRad_showUVH(uvh,cst,pln)
%   matRad_showUVH(uvh,cst,lineStyleIndicator)
%   matRad_showUVH(uvh,cst,pln,lineStyleIndicator)
%
% input
%   uvh:                result struct from fluence optimization/sequencing
%   cst:                matRad cst struct
%   pln:                (now optional) matRad pln struct,
%                       standard uses Dose [Gy]
%   lineStyleIndicator: (optional) integer (1,2,3,4) to indicate the current linestyle
%                       (hint: use different lineStyles to overlay
%                       different uvhs)
%
% output
%   graphical display of UVH   
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

if ~exist('lineStyleIndicator','var') || isempty(lineStyleIndicator)
    lineStyleIndicator = 1;
end

% create new figure and set default line style indicator if not explictly
% specified
hold on;

%reduce cst
visibleIx = cellfun(@(c) c.Visible == 1,cst(:,5));
cstNames = cst(visibleIx,2);
cstInfo = cst(visibleIx,5);
uvh = uvh(visibleIx);

numOfVois = numel(cstNames);
        
%% print the uvh

%try to get colors from cst
try
    colorMx = cellfun(@(c) c.visibleColor,cstInfo,'UniformOutput',false);
    colorMx = cell2mat(colorMx);
catch
    colorMx    = colorcube;
    colorMx    = colorMx(1:floor(64/numOfVois):64,:);
end

lineStyles = {'-',':','--','-.'};

maxUVHvol  = 0;
maxUVHdose = 0;
maxDVHdose = 0;

for i = 1:numOfVois
    % cut off at the first zero value where there is no more signal
    % behind
    currUvhIx = max([1 find(uvh(i).volumePoints>0,1,'last')]);
    currUvh = [uvh(i).doseGrid(1:currUvhIx)*pln.numOfFractions;uvh(i).volumePoints(1:currUvhIx)];
    currDvhIx  = max([1 find(dvh(i).volumePoints>0,1,'last')]);
    currDvh = [dvh(i).doseGrid(1:currDvhIx)*pln.numOfFractions;dvh(i).volumePoints(1:currDvhIx)];
    
    plot(currUvh(1,:),currUvh(2,:),'LineWidth',2,'Color',colorMx(i,:), ...
        'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cstNames{i})
    
    maxUVHvol  = max(maxUVHvol,max(currUvh(2,:)));
    maxUVHdose = max(maxUVHdose,max(currUvh(1,:)));
    maxDVHdose = max(maxDVHdose,max(currDvh(1,:)));
end

fontSizeValue = 14;
myLegend = legend('show','location','NorthEast');
set(myLegend,'FontSize',10,'Interpreter','none');
legend boxoff

if ~exist('doseWindow', 'var') 
    doseWindow = [0 1.4*maxUVHdose];
end

ylim([0 1.05*maxUVHvol]);
xlim(doseWindow);

grid on,grid minor
box(gca,'on');
set(gca,'LineWidth',0.5,'FontSize',fontSizeValue);
ylabel('Volume [%]','FontSize',fontSizeValue)

if exist('pln','var') && ~isempty(pln)
    if strcmp(pln.bioParam.model,'none')
        xlabel('Dose Uncertainty [Gy]','FontSize',fontSizeValue);
    else
        xlabel('RBE x Dose Uncertainty [Gy(RBE)]','FontSize',fontSizeValue);
    end
else
    xlabel('Dose Uncertainty [Gy]','FontSize',fontSizeValue);
end
