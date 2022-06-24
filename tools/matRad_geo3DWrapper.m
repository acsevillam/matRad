function [hCube,hCt] = matRad_geo3DWrapper(axesHandle,ct,cst,cube,window,alpha,contourColorMap,...
                                                                          voiSelection,colorBarLabel,boolPlotLegend,varargin)
% matRad tool function to directly plot a complete slice of a ct with dose
% including contours and isolines.
%
% call
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,thresh)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,alpha)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,contourColorMap)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,doseColorMap)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,doseWindow)
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,doseIsoLevels)
%               ...
% [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSliceWrapper(axesHandle,ct,cst,cubeIdx,dose,plane,slice,thresh,alpha,contourColorMap,...
%                                                                          doseColorMap,doseWindow,doseIsoLevels,voiSelection,colorBarLabel,boolPlotLegend,...)
%
% input (required)
%   axesHandle      handle to axes the slice should be displayed in
%   ct              matRad ct struct
%   cst             matRad cst struct
%   cubeIdx         Index of the desired cube in the ct struct
%   dose            dose cube
%
% input (optional / empty)
%   thresh          threshold for display of dose values
%   alpha           alpha value for the dose overlay
%   contourColorMap colormap for the VOI contours
%   doseColorMap    colormap for the dose
%   doseWindow      dose value window
%   doseIsoLevels   levels defining the isodose contours
%   voiSelection    logicals defining the current selection of contours
%                   that should be plotted. Can be set to [] to plot
%                   all non-ignored contours.
%   colorBarLabel   string defining the yLabel of the colorBar
%   boolPlotLegend  boolean if legend should be plottet or not
%   varargin        additional input parameters that are passed on to
%                   individual plotting functions (e.g. 'LineWidth',1.5)
%   
%
% output
%   hCMap       handle to the colormap
%   hDose       handle to the dose plot
%   hCt         handle to the ct plot
%   hContour    handle to the contour plot
%   hIsoDose    handle to iso dose contours
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

% Handle the argument list

if ~exist('window','var') || isempty(window)
    window = [];
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = [0.1 0.1];
end
if ~exist('contourColorMap','var') || isempty(contourColorMap)
   contourColorMap = [];
end

if ~exist('voiSelection','var') || isempty(voiSelection)
   voiSelection = [];
end

if ~exist('colorBarLabel','var') || isempty(colorBarLabel)
   colorBarLabel = [];
end

if ~exist('boolPlotLegend','var') || isempty(boolPlotLegend)
   boolPlotLegend = false;
end

if ~exist('cst','var') || isempty(cst)
   cst = [];
end

set(axesHandle,'XTick',[],'YTick',[],'Visible','off');
axesHandle.Toolbar.Visible = 'off';

ax1 = axes();
ax1.Toolbar.Visible = 'off';
ax1.Visible = 'off';
ax1.OuterPosition = [0 0 0.9 1];
ax1.View = [60 30];
caxis(ax1,window);
hCube = matRad_plotCube3D(ax1,cube,'jet',window,alpha(1));
axis equal;

ax2 = axes();
ax2.Toolbar.Visible = 'off';
ax2.Color = [0 0 0 1];
ax2.OuterPosition = [0 0 0.9 1];
ax2.View = [60 30];
ax2.Visible = 'on';
xlabel('x [cm]');
ylabel('y [cm] ');
zlabel('z [cm]');
hCt = matRad_plotCube3D(ax2,(ct.cubeHU{1}+1000),'bone',[],alpha(2));
axis equal;

axes(ax1);


% Link two axes together
hlink = linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget'});

% get everthin lined up
cb1 = colorbar(ax1,'Position',[0.86 0.11 0.05 0.815]); % four-elements vector to specify Position [left bottom width height]
cb1.Label.String = colorBarLabel;
cb1.Label.FontSize = 14;

end

