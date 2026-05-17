function matRad_showDVHFromSampling(caSamp, scale, cst, pln, scenarios, doseWindow, dvhType, refScen, lineStyleIndicator, varargin)
% matRad_showDVHFromSampling displays sampled dose-volume histograms
%
% call
%   matRad_showDVHFromSampling(caSamp,scale,cst,pln,scenarios)
%   matRad_showDVHFromSampling(caSamp,scale,cst,pln,scenarios,doseWindow,dvhType,refScen,lineStyleIndicator)
%   matRad_showDVHFromSampling(...,'ctScenProb',ctScenProb)
%   matRad_showDVHFromSampling(...,'scenWeights',scenWeights)
%
% input
%   caSamp:             sampled scenario result struct with per-scenario DVHs
%   scale:              dose scaling factor applied to the DVH dose grid
%   cst:                matRad cst struct
%   pln:                matRad pln struct with multScen scenario model
%   scenarios:          sampled scenario indices to plot
%   doseWindow:         (optional) 1x2 dose axis limits
%   dvhType:            (optional) 'multiscenario', 'minmax', or 'trustband'
%   refScen:            (optional) scenario index highlighted as reference
%   lineStyleIndicator: (optional) integer selecting the line style
%   ctScenProb:         (optional) CT scenario probability override as
%                       Name-Value pair [ctScenId probability]
%   scenWeights:        (optional) resolved scenario weights as Name-Value
%                       pair. If provided, these weights are used directly.
%
% output
%   graphical display of sampled DVH curves or bands
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

matRadCfg = MatRad_Config.instance();
if nargin < 6
    doseWindow = [];
end
if nargin < 7
    dvhType = [];
end
if nargin < 8
    refScen = [];
end
if nargin < 9
    lineStyleIndicator = [];
end

[scale, dvhType, refScen, lineStyleIndicator] = matRad_normalizeSamplingDvhInputs(scale, dvhType, refScen, lineStyleIndicator, matRadCfg);
[ctScenProb, scenWeights] = matRad_parseSamplingDvhOptions(varargin{:});

if ~isempty(ctScenProb) && ~isempty(scenWeights)
    matRadCfg.dispError('Specify either ctScenProb or scenWeights, not both.');
end

plotContext = matRad_buildSamplingDvhPlotContext(cst, scale, refScen, lineStyleIndicator);
hold on;

switch dvhType
    case 'multiscenario'
        [maxDVHvol, maxDVHdose] = matRad_plotMultiScenarioDvh(caSamp, cst, scenarios, plotContext);
    case 'minmax'
        [maxDVHvol, maxDVHdose] = matRad_plotMinMaxSamplingDvh(caSamp, cst, scenarios, plotContext);
    case 'trustband'
        scenWeights = matRad_getSamplingScenarioWeights(pln, numel(caSamp), ctScenProb, scenWeights);
        [maxDVHvol, maxDVHdose] = matRad_plotSamplingDvhTrustband( ...
                                                                  caSamp, cst, scenarios, scenWeights, plotContext);
end

matRad_formatSamplingDvhAxes(pln, doseWindow, maxDVHvol, maxDVHdose);

end

function [scale, dvhType, refScen, lineStyleIndicator] = matRad_normalizeSamplingDvhInputs( ...
                                                                                           scale, dvhType, refScen, lineStyleIndicator, matRadCfg)
if ~exist('dvhType', 'var') || isempty(dvhType)
    dvhType = 'multiscenario';
end

if isstring(dvhType) && isscalar(dvhType)
    dvhType = char(dvhType);
end

validDvhTypes = {'multiscenario', 'minmax', 'trustband'};
if ~ischar(dvhType) || ~any(strcmp(dvhType, validDvhTypes))
    matRadCfg.dispError('Unsupported dvhType. Use ''multiscenario'', ''minmax'', or ''trustband''.');
end

if ~exist('scale', 'var') || isempty(scale)
    scale = 1.0;
end

if ~exist('refScen', 'var') || isempty(refScen)
    refScen = 1;
end

if ~exist('lineStyleIndicator', 'var') || isempty(lineStyleIndicator)
    lineStyleIndicator = 1;
end
end

function [ctScenProb, scenWeights] = matRad_parseSamplingDvhOptions(varargin)
p = inputParser;
p.CaseSensitive = false;
p.addParameter('ctScenProb', [], @(x) isempty(x) || ...
               (isnumeric(x) && ismatrix(x) && size(x, 2) == 2 && ...
                all(isfinite(x(:))) && all(x(:, 2) >= 0)));
p.addParameter('scenWeights', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.parse(varargin{:});
ctScenProb = p.Results.ctScenProb;
scenWeights = p.Results.scenWeights;
end

function plotContext = matRad_buildSamplingDvhPlotContext(cst, scale, refScen, lineStyleIndicator)
cstInfo = cst(:, 5);
numOfVois = size(cst, 1);

try
    colorMx = cellfun(@(c) c.visibleColor, cstInfo, 'UniformOutput', false);
    colorMx = cell2mat(colorMx);
catch
    colorMx = colorcube;
    colorMx = colorMx(1:floor(64 / numOfVois):64, :);
end

plotContext.colorMx = colorMx;
plotContext.lineStyles = {'-', ':', '--', '-.'};
plotContext.scale = scale;
plotContext.refScen = refScen;
plotContext.lineStyle = plotContext.lineStyles{lineStyleIndicator};
end

function [maxDVHvol, maxDVHdose] = matRad_plotMultiScenarioDvh(caSamp, cst, scenarios, plotContext)
maxDVHvol = 0;
maxDVHdose = 0;

for i = 1:size(cst, 1)
    if ~cst{i, 5}.Visible
        continue
    end

    for k = scenarios
        currDvh = matRad_getSamplingDvhCurve(caSamp(k).dvh(i), plotContext.scale);
        pLine = plot(currDvh(1, :), currDvh(2, :), 'LineWidth', 0.5, ...
                     'Color', plotContext.colorMx(i, :), 'LineStyle', plotContext.lineStyle, ...
                     'DisplayName', cst{i, 2});
        matRad_formatReferenceScenarioLine(pLine, k, plotContext.refScen);
        [maxDVHvol, maxDVHdose] = matRad_updateDvhLimits(currDvh, maxDVHvol, maxDVHdose);
    end
end
end

function [maxDVHvol, maxDVHdose] = matRad_plotMinMaxSamplingDvh(caSamp, cst, scenarios, plotContext)
maxDVHvol = 0;
maxDVHdose = 0;

for i = 1:size(cst, 1)
    if ~cst{i, 5}.Visible
        continue
    end

    [allDvh, doseGrid] = matRad_collectSamplingDvhValues(caSamp, i, scenarios, plotContext.scale);
    matRad_plotReferenceSamplingDvh(caSamp, cst, i, plotContext);
    bandDvh = [doseGrid; min(allDvh); max(allDvh)];
    matRad_plotSamplingDvhBand(bandDvh, plotContext.colorMx(i, :));
    [maxDVHvol, maxDVHdose] = matRad_updateDvhLimits(bandDvh([1 3], :), maxDVHvol, maxDVHdose);
end
end

function [maxDVHvol, maxDVHdose] = matRad_plotSamplingDvhTrustband(caSamp, cst, scenarios, scenWeights, plotContext)
maxDVHvol = 0;
maxDVHdose = 0;

for i = 1:size(cst, 1)
    if ~cst{i, 5}.Visible
        continue
    end

    [allDvh, doseGrid] = matRad_collectSamplingDvhValues(caSamp, i, scenarios, plotContext.scale);
    currScenWeights = matRad_getCurrentScenarioWeights(scenWeights, scenarios, numel(caSamp));
    meanDVHVolW = matRad_weightedMean(allDvh, currScenWeights);
    stdDVHVolW = matRad_weightedStd(allDvh, currScenWeights);
    bandDvh = [doseGrid; max(meanDVHVolW - stdDVHVolW, 0); min(meanDVHVolW + stdDVHVolW, 100)];

    matRad_plotSamplingDvhBand(bandDvh, plotContext.colorMx(i, :));
    matRad_plotReferenceSamplingDvh(caSamp, cst, i, plotContext);
    pLine = plot(doseGrid, meanDVHVolW, 'LineWidth', 2, ...
                 'Color', plotContext.colorMx(i, :), 'LineStyle', plotContext.lineStyles{2}, ...
                 'DisplayName', cst{i, 2});
    pLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
    [maxDVHvol, maxDVHdose] = matRad_updateDvhLimits(bandDvh([1 3], :), maxDVHvol, maxDVHdose);
end
end

function [allDvh, doseGrid] = matRad_collectSamplingDvhValues(caSamp, voiIx, scenarios, scale)
doseGrid = caSamp(1).dvh(voiIx).doseGrid * scale;
allDvh = zeros(numel(scenarios), numel(doseGrid));

for s = 1:numel(scenarios)
    allDvh(s, :) = caSamp(scenarios(s)).dvh(voiIx).volumePoints;
end
end

function matRad_plotReferenceSamplingDvh(caSamp, cst, voiIx, plotContext)
refScen = plotContext.refScen;
currDvh = matRad_getSamplingDvhCurve(caSamp(refScen).dvh(voiIx), plotContext.scale);
pLine = plot(currDvh(1, :), currDvh(2, :), 'LineWidth', 2, ...
             'Color', plotContext.colorMx(voiIx, :), 'LineStyle', plotContext.lineStyle, ...
             'DisplayName', cst{voiIx, 2});
pLine.Annotation.LegendInformation.IconDisplayStyle = 'on';
end

function currDvh = matRad_getSamplingDvhCurve(dvh, scale)
ix = matRad_lastDvhPlotIndex(dvh.volumePoints);
currDvh = [dvh.doseGrid(1:ix) * scale; dvh.volumePoints(1:ix)];
end

function matRad_formatReferenceScenarioLine(pLine, scenarioIx, refScen)
if scenarioIx == refScen
    pLine.LineWidth = 2;
    pLine.Annotation.LegendInformation.IconDisplayStyle = 'on';
else
    pLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
end

function matRad_plotSamplingDvhBand(bandDvh, color)
pFill = fill([bandDvh(1, :) flip(bandDvh(1, :))], ...
             [bandDvh(2, :) flip(bandDvh(3, :))], color);
pFill.FaceAlpha = 0.2;
pFill.FaceColor = color;
pFill.EdgeColor = color;
pFill.LineWidth = 0.05;
pFill.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

function [maxDVHvol, maxDVHdose] = matRad_updateDvhLimits(currDvh, maxDVHvol, maxDVHdose)
maxDVHvol = max(maxDVHvol, max(currDvh(2, :)));
maxDVHdose = max(maxDVHdose, max(currDvh(1, :)));
end

function matRad_formatSamplingDvhAxes(pln, doseWindow, maxDVHvol, maxDVHdose)
fontSizeValue = 14;
myLegend = legend('show', 'location', 'NorthEast');
set(myLegend, 'FontSize', 10, 'Interpreter', 'none');
legend boxoff;

doseWindow = matRad_getSamplingDvhDoseWindow(doseWindow, maxDVHdose);
if maxDVHvol <= 0
    maxDVHvol = 1;
end

ylim([0 1.05 * maxDVHvol]);
xlim(doseWindow(1:2));
grid on;
grid minor;
box(gca, 'on');
set(gca, 'LineWidth', .5, 'FontSize', fontSizeValue);
ylabel('Volume [%]', 'FontSize', fontSizeValue);

if strcmp(matRad_getPlanBioModelName(pln), 'none')
    xlabel('Dose [Gy]', 'FontSize', fontSizeValue);
else
    xlabel('RBE x Dose [Gy(RBE)]', 'FontSize', fontSizeValue);
end
end

function doseWindow = matRad_getSamplingDvhDoseWindow(doseWindow, maxDVHdose)
if ~exist('doseWindow', 'var') || isempty(doseWindow)
    doseWindow = [0 1.4 * maxDVHdose];
else
    doseWindow = doseWindow(:)';
end

if numel(doseWindow) >= 2 && all(isfinite(doseWindow(1:2))) && doseWindow(2) > doseWindow(1)
    doseWindow = doseWindow(1:2);
    return
end

fallbackMaxDose = 1.4 * maxDVHdose;
if ~isfinite(fallbackMaxDose) || fallbackMaxDose <= 0
    fallbackMaxDose = 1;
end
doseWindow = [0 fallbackMaxDose];
end

function currScenWeights = matRad_getCurrentScenarioWeights(scenWeights, scenarios, numScenarios)
matRad_cfg = MatRad_Config.instance();

if numel(scenWeights) == numScenarios
    currScenWeights = scenWeights(scenarios(:));
elseif numel(scenWeights) == numel(scenarios)
    currScenWeights = scenWeights(:);
else
    matRad_cfg.dispError('Number of scenario weights does not match selected scenarios.');
end
end

function ix = matRad_lastDvhPlotIndex(volumePoints)
lastPositiveIx = find(volumePoints > 0, 1, 'last');
if isempty(lastPositiveIx)
    ix = 1;
else
    ix = min(lastPositiveIx + 1, numel(volumePoints));
end
end

function meanValue = matRad_weightedMean(values, weights)
weights = weights(:);
meanValue = weights' * values ./ sum(weights);
end

function stdValue = matRad_weightedStd(values, weights)
meanValue = matRad_weightedMean(values, weights);
stdValue = sqrt(matRad_weightedMean(bsxfun(@minus, values, meanValue).^2, weights));
end

function modelName = matRad_getPlanBioModelName(pln)
modelName = 'none';

if isfield(pln, 'bioModel')
    bioModel = pln.bioModel;
elseif isfield(pln, 'bioParam')
    bioModel = pln.bioParam;
else
    return
end

if ischar(bioModel)
    modelName = bioModel;
elseif isstring(bioModel) && isscalar(bioModel)
    modelName = char(bioModel);
elseif isstruct(bioModel) && isfield(bioModel, 'model')
    modelName = bioModel.model;
elseif isobject(bioModel) && isprop(bioModel, 'model')
    modelName = bioModel.model;
end
end
