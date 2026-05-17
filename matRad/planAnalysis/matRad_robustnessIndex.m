function [robustnessCube, robPassRate, robustnessFig] = matRad_robustnessIndex(meanCube, stdCube, refDose, criteria, method, ct, cst, slice)
% matRad_robustnessIndex calculates a target robustness index
%
% call
%   [robustnessCube,robPassRate] = matRad_robustnessIndex(meanCube,stdCube,refDose,criteria,method,ct,cst)
%   [robustnessCube,robPassRate,robustnessFig] = matRad_robustnessIndex(meanCube,stdCube,refDose,criteria,method,ct,cst,slice)
%
% input
%   meanCube:         expected dose cube in per-fraction units
%   stdCube:          dose standard deviation cube in per-fraction units
%   refDose:          per-fraction reference dose
%   criteria:         1x2 vector [mean dose deviation %, std dose %]
%   method:           robustness index method, 'index1' or 'index2'
%   ct:               matRad ct struct
%   cst:              matRad cst cell array used to identify TARGET voxels
%   slice:            optional CT slice for visualization
%
% output
%   robustnessCube:   robustness index cube
%   robPassRate:      target voxel pass rate in percent
%   robustnessFig:    optional robustness figure
%
% References
%   [1] Sevilla-Moreno AC, Puerta-Yepes ME, Wahl N, Benito-Herce R,
%       Cabal-Arango G. Interval Analysis-Based Optimization: A Robust Model
%       for Intensity-Modulated Radiotherapy (IMRT). Cancers (Basel).
%       2025 Feb 3;17(3):504. doi:10.3390/cancers17030504.
%   [2] Cabal G, Jaekel O. Dynamic Target Definition: A novel approach for
%       PTV definition in ion beam therapy. Radiother Oncol. 2013;107:227-233.
%       doi:10.1016/j.radonc.2013.03.010.
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

if nargin < 8
    slice = [];
end

matRad_validateRobustnessIndexInput(meanCube, stdCube, refDose, criteria);
method = matRad_normalizeRobustnessIndexMethod(method);

targetMask = matRad_getRobustnessTargetMask(meanCube, ct, cst);
meanDoseThreshold = criteria(1);
stdThreshold = criteria(2);

meanDoseDeviationPct = abs(meanCube - refDose) ./ refDose * 100;
stdPct = stdCube ./ refDose * 100;

switch method
    case 'index1'
        robustnessCube = sqrt((meanDoseDeviationPct ./ meanDoseThreshold).^2 + ...
                              (stdPct ./ stdThreshold).^2);
        passMask = robustnessCube < 1;
    case 'index2'
        meanDoseCrit = meanDoseDeviationPct <= meanDoseThreshold;
        stdCrit = stdPct <= stdThreshold;
        robustnessCube = meanDoseCrit & stdCrit;
        passMask = robustnessCube == 1;
end

robPassRate = matRad_calcRobustnessPassRate(passMask, targetMask);
robustnessFig = matRad_plotRobustnessIndex(robustnessCube, targetMask, robPassRate, ...
                                           criteria, method, ct, cst, slice);

end

function matRad_validateRobustnessIndexInput(meanCube, stdCube, refDose, criteria)
validInput = isequal(size(meanCube), size(stdCube)) && ...
    isnumeric(refDose) && isscalar(refDose) && isfinite(refDose) && refDose > 0 && ...
    isnumeric(criteria) && numel(criteria) == 2 && all(isfinite(criteria)) && ...
    all(criteria > 0);

if ~validInput
    matRadCfg = MatRad_Config.instance();
    matRadCfg.dispError('Invalid robustness index input.');
end
end

function method = matRad_normalizeRobustnessIndexMethod(method)
if isnumeric(method) && isscalar(method)
    if method == 1
        method = 'index1';
    elseif method == 2
        method = 'index2';
    end
end

if isstring(method) && isscalar(method)
    method = char(method);
end

if ischar(method)
    method = lower(method);
end

if ~(ischar(method) && any(strcmp(method, {'index1', 'index2'})))
    matRadCfg = MatRad_Config.instance();
    matRadCfg.dispError('Robustness index method must be ''index1'' or ''index2''.');
end
end

function targetMask = matRad_getRobustnessTargetMask(cube, ct, cst)
matRadCfg = MatRad_Config.instance();

if isempty(cst) || ~iscell(cst)
    matRadCfg.dispError('Robustness index calculation requires a cst input.');
end

refScen = matRad_getRobustnessReferenceScenario(ct);
targetRows = find(strcmpi(cst(:, 3), 'TARGET'))';
targetMask = false(size(cube));

for i = targetRows
    if numel(cst{i, 4}) >= refScen && ~isempty(cst{i, 4}{refScen})
        targetMask(cst{i, 4}{refScen}) = true;
    end
end

if ~any(targetMask(:))
    matRadCfg.dispError('Robustness index calculation requires at least one TARGET voxel.');
end
end

function refScen = matRad_getRobustnessReferenceScenario(ct)
if isstruct(ct) && isfield(ct, 'refScen')
    refScen = ct.refScen;
else
    refScen = 1;
end
end

function robPassRate = matRad_calcRobustnessPassRate(passMask, targetMask)
numTargetVoxels = nnz(targetMask);
if numTargetVoxels > 0
    robPassRate = 100 * nnz(passMask(targetMask)) / numTargetVoxels;
else
    robPassRate = NaN;
end
end

function robustnessFig = matRad_plotRobustnessIndex(robustnessCube, targetMask, robPassRate, criteria, method, ct, cst, slice)
robustnessFig = [];
if isempty(slice)
    return
end

plotCube = double(robustnessCube);
plotCube(~targetMask) = NaN;

robustnessFig = figure;
set(robustnessFig, 'Color', [1 1 1]);
if matRad_canUseRobustnessSliceWrapper(ct)
    refScen = matRad_getRobustnessReferenceScenario(ct);
    matRad_plotSliceWrapper(gca, ct, cst, refScen, plotCube, 3, slice, [], [], ...
                            colorcube, matRad_getRobustnessColormap(method), ...
                            matRad_getRobustnessDoseWindow(method), [], [], ...
                            matRad_getRobustnessColorbarLabel(method), [], ...
                            'LineWidth', 1.5);
else
    imagesc(plotCube(:, :, slice), matRad_getRobustnessDoseWindow(method));
    colormap(gca, matRad_getRobustnessColormap(method));
    colorbar;
end

title(sprintf('RI = %.4f (%g%% / %g%%)', robPassRate / 100, criteria(1), criteria(2)));
matRad_addRobustnessMetricTextbox(gca, robPassRate);
end

function tf = matRad_canUseRobustnessSliceWrapper(ct)
tf = isstruct(ct) && isfield(ct, 'cubeHU') && iscell(ct.cubeHU) && ...
    isfield(ct, 'resolution') && (isfield(ct, 'cubeDim') || isfield(ct, 'dimensions'));
end

function colorMap = matRad_getRobustnessColormap(method)
switch method
    case 'index1'
        colorMap = matRad_getRobustnessIndex1ColorMap();
    otherwise
        colorMap = [];
end
end

function colorMap = matRad_getRobustnessIndex1ColorMap()
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

function doseWindow = matRad_getRobustnessDoseWindow(method)
switch method
    case 'index1'
        doseWindow = [0 5.01];
    otherwise
        doseWindow = [0 2];
end
end

function label = matRad_getRobustnessColorbarLabel(method)
switch method
    case 'index1'
        label = 'Delta Index';
    otherwise
        label = [];
end
end

function matRad_addRobustnessMetricTextbox(axesHandle, robPassRate)
delete(findobj(axesHandle, 'Tag', 'matRadMetricTextbox'));
text(axesHandle, 0.02, 0.04, matRad_formatRobustnessIndexText(robPassRate), ...
     'Units', 'normalized', ...
     'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'bottom', ...
     'BackgroundColor', [1 1 1], ...
     'EdgeColor', [0.2 0.2 0.2], ...
     'Margin', 4, ...
     'FontWeight', 'bold', ...
     'Tag', 'matRadMetricTextbox');
end

function metricText = matRad_formatRobustnessIndexText(robPassRate)
if isempty(robPassRate) || ~isnumeric(robPassRate) || ~isscalar(robPassRate) || ...
        ~isfinite(robPassRate)
    metricText = 'RI = n/a';
else
    metricText = sprintf('RI = %.4f', robPassRate / 100);
end
end
