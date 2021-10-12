function matRad_showDVH_sampled(dvh,cst,pln,doseWindow,lineStyleIndicator)
% matRad dvh visualizaion
% 
% call
%   matRad_showDVH(dvh,cst,pln,lineStyleIndicator)
%
% input
%   result:             result struct from fluence optimization/sequencing
%   cst:                matRad cst struct
%   pln:                matRad pln struct
%   lineStyleIndicator: integer (1,2,3,4) to indicate the current linestyle
%                       (hint: use different lineStyles to overlay
%                       different dvhs)
%
% output
%   graphical display of DVH   
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

numOfVois = numel(cstNames);
        
%% print the dvh
try
    colorMx = cellfun(@(c) c.visibleColor,cstInfo,'UniformOutput',false);
    colorMx = cell2mat(colorMx);
catch
    colorMx    = colorcube;
    colorMx    = colorMx(1:floor(64/numOfVois):64,:);
end

lineStyles = {'-',':','--','-.'};

maxDVHvol  = 0;
maxDVHdose = 0;

for i = 1:numOfVois
    if cst{i,5}.Visible
        % cut off at the first zero value where there is no more signal
        % behind
        currIx      = max([1 find(dvh{1,1}(i).volumePoints>0,1,'last')]);
        currDvh = [dvh{1,1}(i).doseGrid(1:currIx)*pln.numOfFractions;dvh{1,1}(i).volumePoints(1:currIx)];
        lowerIx      = max([1 find(dvh{1,2}(i).volumePoints>0,1,'last')]);
        lowerDvh = [dvh{1,2}(i).doseGrid(1:lowerIx)*pln.numOfFractions;dvh{1,2}(i).volumePoints(1:lowerIx)];
        upperIx      = max([1 find(dvh{1,3}(i).volumePoints>0,1,'last')]);
        upperDvh = [dvh{1,3}(i).doseGrid(1:upperIx)*pln.numOfFractions;dvh{1,3}(i).volumePoints(1:upperIx)];

        plot(currDvh(1,:),currDvh(2,:),'LineWidth',2,'Color',colorMx(i,:), ...
            'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cst{i,2})

        f = fill([lowerDvh(1,:) flip(upperDvh(1,:))],[lowerDvh(2,:) flip(upperDvh(2,:))],colorMx(i,:));
        f.FaceAlpha = 0.2;
        f.FaceColor = colorMx(i,:);
        f.EdgeColor = colorMx(i,:);
        f.LineWidth = 0.05;
        f.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        
        maxDVHvol  = max(maxDVHvol,max(currDvh(2,:)));
        maxDVHdose = max(maxDVHdose,max(currDvh(1,:)));
    end
end

fontSizeValue = 14;
myLegend = legend('show','location','NorthEast');
set(myLegend,'FontSize',10,'Interpreter','none');
legend boxoff

if ~exist('doseWindow', 'var') 
    doseWindow = [0 1.4*maxDVHdose];
end

ylim([0 1.05*maxDVHvol]);
xlim(doseWindow);

grid on,grid minor
box(gca,'on');
set(gca,'LineWidth',.5,'FontSize',fontSizeValue);
ylabel('Volume [%]','FontSize',fontSizeValue)

if strcmp(pln.bioParam.model,'none')
     xlabel('Dose [Gy]','FontSize',fontSizeValue);
else
     xlabel('RBE x Dose [Gy(RBE)]','FontSize',fontSizeValue);
end