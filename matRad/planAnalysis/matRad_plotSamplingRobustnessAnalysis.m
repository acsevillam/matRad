function robustnessFig = matRad_plotSamplingRobustnessAnalysis(robustnessAnalysis, ct, cst, slice, varargin)
% matRad_plotSamplingRobustnessAnalysis plots a sampling robustness index
%
% call
%   robustnessFig = matRad_plotSamplingRobustnessAnalysis(robustnessAnalysis,ct,cst,slice)
%
% input
%   robustnessAnalysis: sampling robustness analysis struct
%   ct:                 matRad ct struct
%   cst:                matRad cst cell array
%   slice:              CT slice to plot
%   varargin:           optional Name/Value pairs:
%                       - 'axesHandle': axes used for plotting
%                       - 'method': 'index1' or 'index2'
%                       - 'plane': matRad plane index
%                       - 'contourColorMap': VOI contour colors
%
% output
%   robustnessFig:      figure handle when a new figure is created
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4 || isempty(slice)
    robustnessFig = [];
    return
end

p = inputParser;
p.CaseSensitive = false;
p.addParameter('axesHandle', [], @(axesHandle) isempty(axesHandle) || ...
               strcmp(get(axesHandle, 'type'), 'axes'));
p.addParameter('method', 'index1', @(method) ischar(method) || ...
               (isstring(method) && isscalar(method)));
p.addParameter('plane', 3, @(plane) isnumeric(plane) && isscalar(plane) && any(plane == [1 2 3]));
p.addParameter('contourColorMap', [], @(colorMap) isempty(colorMap) || ...
               (isnumeric(colorMap) && size(colorMap, 2) == 3));
parse(p, varargin{:});

method = matRad_normalizeSamplingRobustnessMethod(p.Results.method);
axesHandle = p.Results.axesHandle;
if isempty(axesHandle)
    robustnessFig = figure('Name', 'Sampling robustness analysis');
    set(robustnessFig, 'Color', [1 1 1]);
    axesHandle = gca;
else
    robustnessFig = ancestor(axesHandle, 'figure');
end

[plotCube, plotWindow, colorMap, colorBarLabel] = ...
    matRad_getSamplingRobustnessPlotData(robustnessAnalysis, method);

matRad_plotSlice(ct, 'axesHandle', axesHandle, 'cst', cst, 'cubeIdx', 1, ...
                 'dose', plotCube, 'plane', p.Results.plane, 'slice', slice, ...
                 'contourColorMap', p.Results.contourColorMap, ...
                 'doseColorMap', colorMap, 'doseWindow', plotWindow, ...
                 'colorBarLabel', colorBarLabel);

title(axesHandle, matRad_getSamplingRobustnessPlotTitle(robustnessAnalysis, method));
matRad_addSamplingRobustnessMetricTextbox(axesHandle, robustnessAnalysis, method);

end

function method = matRad_normalizeSamplingRobustnessMethod(method)
method = lower(char(method));
if ~any(strcmp(method, {'index1', 'index2'}))
    matRadCfg = MatRad_Config.instance();
    matRadCfg.dispError('Sampling robustness plot method must be ''index1'' or ''index2''.');
end
end

function [plotCube, plotWindow, colorMap, colorBarLabel] = ...
    matRad_getSamplingRobustnessPlotData(robustnessAnalysis, method)
switch method
    case 'index1'
        plotCube = robustnessAnalysis.index1.robustnessCube;
        plotWindow = [0 5.01];
        colorMap = matRad_getSamplingRobustnessIndex1ColorMap();
        colorBarLabel = 'Delta Index';
    case 'index2'
        plotCube = double(robustnessAnalysis.index2.robustnessCube);
        plotWindow = [0 2.01];
        colorMap = [];
        colorBarLabel = [];
end

if isfield(robustnessAnalysis, 'evaluableTargetMask') && ...
        isequal(size(plotCube), size(robustnessAnalysis.evaluableTargetMask))
    targetMask = robustnessAnalysis.evaluableTargetMask;
    plotCube(~targetMask) = NaN;
end
end

function titleText = matRad_getSamplingRobustnessPlotTitle(robustnessAnalysis, method)
robustnessIndex = robustnessAnalysis.(method).robustnessIndex;
if isempty(robustnessIndex) || ~isfinite(robustnessIndex)
    titleText = sprintf('%s unavailable', method);
else
    titleText = sprintf('%s RI %.3f', method, robustnessIndex);
end
end

function matRad_addSamplingRobustnessMetricTextbox(axesHandle, robustnessAnalysis, method)
delete(findobj(axesHandle, 'Tag', 'matRadMetricTextbox'));
text(axesHandle, 0.02, 0.04, matRad_formatSamplingRobustnessIndexText( ...
                                                                      robustnessAnalysis.(method).robustnessIndex), ...
     'Units', 'normalized', ...
     'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'bottom', ...
     'BackgroundColor', [1 1 1], ...
     'EdgeColor', [0.2 0.2 0.2], ...
     'Margin', 4, ...
     'FontWeight', 'bold', ...
     'Tag', 'matRadMetricTextbox');
end

function metricText = matRad_formatSamplingRobustnessIndexText(robustnessIndex)
if isempty(robustnessIndex) || ~isnumeric(robustnessIndex) || ...
        ~isscalar(robustnessIndex) || ~isfinite(robustnessIndex)
    metricText = 'RI = n/a';
else
    metricText = sprintf('RI = %.4f', robustnessIndex);
end
end

function colorMap = matRad_getSamplingRobustnessIndex1ColorMap()
numColors = 256;
maxRobustnessIndex = 5.01;
numPassColors = round(1 / maxRobustnessIndex * numColors);
numFailColors = numColors - numPassColors;

passColorMap = [linspace(0.40, 1, numPassColors)' ...
                ones(numPassColors, 1) ...
                linspace(0.40, 1, numPassColors)'];
failColorMap = matRad_getColormap('gammaIndex', 2 * numFailColors);
colorMap = [passColorMap; failColorMap(numFailColors + 1:end - 1, :)];
end
