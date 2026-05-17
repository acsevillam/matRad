function expectedDoseDifferenceFig = matRad_plotExpectedDoseDifferenceAnalysis(expectedDoseDifferenceAnalysis, ct, cst, slice, varargin)
% matRad_plotExpectedDoseDifferenceAnalysis plots signed expected dose difference
%
% call
%   expectedDoseDifferenceFig = matRad_plotExpectedDoseDifferenceAnalysis(expectedDoseDifferenceAnalysis,ct,cst,slice)
%
% input
%   expectedDoseDifferenceAnalysis: expected dose difference analysis struct
%   ct:                  matRad ct struct
%   cst:                 matRad cst cell array
%   slice:               CT slice to plot
%   varargin:            optional Name/Value pairs:
%                        - 'axesHandle': axes used for plotting
%                        - 'plane': matRad plane index
%                        - 'doseWindow': signed dose-difference color window
%                        - 'contourColorMap': VOI contour colors
%
% output
%   expectedDoseDifferenceFig:      figure handle when a new figure is created
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
    expectedDoseDifferenceFig = [];
    return
end

p = inputParser;
p.CaseSensitive = false;
p.addParameter('axesHandle', [], @(axesHandle) isempty(axesHandle) || ...
               strcmp(get(axesHandle, 'type'), 'axes'));
p.addParameter('plane', 3, @(plane) isnumeric(plane) && isscalar(plane) && any(plane == [1 2 3]));
p.addParameter('doseWindow', [], @(window) isempty(window) || ...
               (isnumeric(window) && numel(window) == 2 && window(1) < window(2)));
p.addParameter('contourColorMap', [], @(colorMap) isempty(colorMap) || ...
               (isnumeric(colorMap) && size(colorMap, 2) == 3));
parse(p, varargin{:});

axesHandle = p.Results.axesHandle;
if isempty(axesHandle)
    expectedDoseDifferenceFig = figure('Name', 'Expected dose difference analysis');
    set(expectedDoseDifferenceFig, 'Color', [1 1 1]);
    axesHandle = axes('Parent', expectedDoseDifferenceFig);
else
    expectedDoseDifferenceFig = ancestor(axesHandle, 'figure');
end

[plotCube, plotWindow, colorBarLabel] = ...
    matRad_getExpectedDoseDifferencePlotData(expectedDoseDifferenceAnalysis, p.Results.doseWindow);
plotCube = matRad_reshapeExpectedDoseDifferenceCubeForPlot(plotCube, ct);
expectedDoseDifferenceColorMap = matRad_getExpectedDoseDifferenceColorMap(128);
if matRad_canUsePlotSliceForExpectedDoseDifference(ct, plotCube)
    matRad_plotSlice(ct, 'axesHandle', axesHandle, 'cst', cst, 'cubeIdx', 1, ...
                     'dose', plotCube, 'plane', p.Results.plane, 'slice', slice, ...
                     'contourColorMap', p.Results.contourColorMap, ...
                     'doseColorMap', expectedDoseDifferenceColorMap, ...
                     'doseWindow', plotWindow, ...
                     'colorBarLabel', colorBarLabel);
else
    imagesc(axesHandle, matRad_getExpectedDoseDifferenceSlice(plotCube, p.Results.plane, slice), ...
            plotWindow);
    axis(axesHandle, 'image');
    colormap(axesHandle, expectedDoseDifferenceColorMap);
    hColorBar = colorbar(axesHandle);
    set(get(hColorBar, 'YLabel'), 'String', colorBarLabel);
end

title(axesHandle, sprintf('%s: %s', expectedDoseDifferenceAnalysis.referenceName, ...
                          expectedDoseDifferenceAnalysis.status));

end

function [plotCube, plotWindow, colorBarLabel] = ...
    matRad_getExpectedDoseDifferencePlotData(expectedDoseDifferenceAnalysis, doseWindow)
if isfield(expectedDoseDifferenceAnalysis, 'signedExpectedDoseDifferenceCube')
    plotCube = expectedDoseDifferenceAnalysis.signedExpectedDoseDifferenceCube;
    colorBarLabel = ['E[D - ref] ' ...
                     matRad_getExpectedDoseDifferenceUnitLabel(expectedDoseDifferenceAnalysis)];
else
    plotCube = expectedDoseDifferenceAnalysis.signedReferenceProbabilityCube;
    colorBarLabel = 'P(D > ref) - P(D < ref)';
end

if isempty(doseWindow) && isfield(expectedDoseDifferenceAnalysis, 'doseWindow')
    plotWindow = expectedDoseDifferenceAnalysis.doseWindow;
elseif isempty(doseWindow)
    plotWindow = matRad_getSymmetricFiniteWindow(plotCube);
else
    plotWindow = doseWindow;
end
end

function plotWindow = matRad_getSymmetricFiniteWindow(plotCube)
finiteValues = plotCube(isfinite(plotCube));
if isempty(finiteValues)
    plotWindow = [-1 1];
    return
end

windowAbs = max(abs(finiteValues(:)));
if windowAbs == 0
    windowAbs = 1;
end
plotWindow = [-windowAbs windowAbs];
end

function colorMap = matRad_getExpectedDoseDifferenceColorMap(numColors)
if nargin < 1 || isempty(numColors)
    numColors = 128;
end

numLower = floor(numColors / 2);
numUpper = numColors - numLower;

blueToWhite = [linspace(0, 1, numLower)' ...
               linspace(0, 1, numLower)' ...
               ones(numLower, 1)];
whiteToRed = [ones(numUpper, 1) ...
              linspace(1, 0, numUpper)' ...
              linspace(1, 0, numUpper)'];
colorMap = [blueToWhite; whiteToRed];
end

function tf = matRad_canUsePlotSliceForExpectedDoseDifference(ct, plotCube)
tf = isfield(ct, 'cubeDim') && numel(size(plotCube)) == numel(ct.cubeDim) && ...
    all(size(plotCube) == ct.cubeDim);
end

function sliceData = matRad_getExpectedDoseDifferenceSlice(plotCube, plane, slice)
if ndims(plotCube) < 3
    sliceData = plotCube;
    return
end

switch plane
    case 1
        sliceData = squeeze(plotCube(slice, :, :));
    case 2
        sliceData = squeeze(plotCube(:, slice, :));
    otherwise
        sliceData = squeeze(plotCube(:, :, slice));
end
end

function unitLabel = matRad_getExpectedDoseDifferenceUnitLabel(expectedDoseDifferenceAnalysis)
quantity = '';
if isfield(expectedDoseDifferenceAnalysis, 'quantity')
    quantity = char(expectedDoseDifferenceAnalysis.quantity);
end

if strncmp(quantity, 'RBExD', 5) || strncmp(quantity, 'RBExDose', 8)
    unitLabel = '[Gy(RBE)]';
else
    unitLabel = '[Gy]';
end
end

function plotCube = matRad_reshapeExpectedDoseDifferenceCubeForPlot(plotCube, ct)
if isfield(ct, 'cubeDim') && numel(plotCube) == prod(ct.cubeDim) && ...
        ~isequal(size(plotCube), ct.cubeDim)
    plotCube = reshape(plotCube, ct.cubeDim);
end
end
