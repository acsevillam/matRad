function [robustnessCube,robPassRate,robustnessFig] = matRad_robustnessIndex(meanCube,stdCube,refDose,criteria,method,ct,cst,slice,varargin)
% matRad_robustnessIndex calculates a sampling robustness index
%
% call
%   [robustnessCube,robPassRate] = matRad_robustnessIndex(meanCube,stdCube,refDose,criteria,method,ct,cst)
%   [robustnessCube,robPassRate,robustnessFig] = matRad_robustnessIndex(meanCube,stdCube,refDose,criteria,method,ct,cst,slice)
%   robustnessFig = matRad_robustnessIndex('plotAggregate',robustnessCube,targetMask,robPassRate,method,criteria,slice,ct,cst,'plotMode',plotMode)
%
% input
%   meanCube:          expected dose cube
%   stdCube:           dose standard deviation cube
%   refDose:           reference dose used for percentage criteria
%   criteria:          1x2 vector [mean dose deviation %, std dose %]
%   method:            robustness index method, 'index1' for the combined
%                      method in [1] or 'index2' for the binary criterion
%                      derived from the tCTV(m,v) definition in [2]
%   ct:                (optional) matRad ct struct
%   cst:               matRad cst struct used to identify TARGET voxels
%   slice:             (optional) CT slice used to create a figure
%   varargin:          optional cst input for contours in 'plotAggregate'
%                      mode. The aggregate target mask is supplied directly
%                      in that mode. Optional 'plotMode' for index1 can be
%                      'continuous' or 'binary'. Default is 'continuous'.
%                      Optional 'robustnessTargetMode' can be 'all',
%                      'include', or 'exclude'. Optional
%                      'robustnessTargets' accepts target names or cst row
%                      indices.
%
% output
%   robustnessCube:    robustness cube for the selected method
%   robPassRate:       percentage of target voxels passing the criterion
%   robustnessFig:     robustness figure if slice was provided
%
% methods
%   index1: computes the Delta Index, a combined normalized mean-dose
%           deviation and standard-deviation robustness metric. Target
%           voxels pass when
%           sqrt((meanDoseDeviationPct/meanDoseThreshold)^2 +
%           (stdPct/stdThreshold)^2) < 1 [1].
%   index2: computes a binary criterion. Target voxels pass when both
%           meanDoseDeviationPct <= meanDoseThreshold and
%           stdPct <= stdThreshold [2].
%
% References
%   [1] Sevilla-Moreno AC, Puerta-Yepes ME, Wahl N, Benito-Herce R,
%       Cabal-Arango G. Interval Analysis-Based Optimization: A Robust Model
%       for Intensity-Modulated Radiotherapy (IMRT). Cancers (Basel).
%       2025 Feb 3;17(3):504. doi:10.3390/cancers17030504.
%       PMID:39941871; PMCID:PMC11816179.
%       https://pmc.ncbi.nlm.nih.gov/articles/PMC11816179/
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

if nargin < 6
    ct = [];
end

if nargin < 7
    cst = [];
end

if nargin < 8
    slice = [];
end

if isRobustnessIndexAction(meanCube)
    [robustnessCube,robPassRate,robustnessFig] = runRobustnessIndexAction( ...
        meanCube,stdCube,refDose,criteria,method,ct,cst,slice,varargin{:});
    return;
end

validateRobustnessInput(meanCube,stdCube,refDose,criteria);
method = normalizeRobustnessMethod(method);

meanDoseThreshold = criteria(1);
stdThreshold = criteria(2);

robustnessOptions = parseRobustnessOptions(varargin);
[targetMask,doseMask,refScen] = getRobustnessMasks(meanCube,stdCube,ct,cst, ...
    robustnessOptions.robustnessTargetMode,robustnessOptions.robustnessTargets);
plotMode = robustnessOptions.plotMode;

meanDoseDeviationPct = abs(meanCube - refDose)./refDose*100;
stdPct = stdCube./refDose*100;

switch method
    case 'index1'
        robustnessCube = sqrt((meanDoseDeviationPct./meanDoseThreshold).^2 + ...
            (stdPct./stdThreshold).^2);
        robPassRate = calcPassRate(robustnessCube < 1,targetMask);
        robustnessFig = plotRobustnessIndex1(robustnessCube,targetMask, ...
            robPassRate,meanDoseThreshold,stdThreshold,slice,ct,cst,refScen,plotMode);
    case 'index2'
        meanDoseCrit = meanDoseDeviationPct <= meanDoseThreshold;
        stdCrit = stdPct <= stdThreshold;
        robustnessCube = meanDoseCrit & stdCrit;
        robPassRate = calcPassRate(robustnessCube,targetMask);
        robustnessFig = plotRobustnessIndex2(robustnessCube,targetMask,doseMask, ...
            robPassRate,meanDoseThreshold,stdThreshold,slice,ct,cst,refScen);
end

end

function tf = isRobustnessIndexAction(value)
tf = (ischar(value) || (isstring(value) && isscalar(value))) && ...
    any(strcmpi(char(value),{'plotAggregate'}));
end

function [robustnessCube,robPassRate,robustnessFig] = runRobustnessIndexAction(action,robustnessCubeIn,targetMask,robPassRateIn,method,criteria,slice,ct,varargin)
switch lower(char(action))
    case 'plotaggregate'
        plotArgs = varargin;
        if isempty(varargin)
            cst = [];
        elseif isRobustnessPlotMode(varargin{1}) || isPlotModeName(varargin{1})
            cst = [];
        else
            cst = varargin{1};
            plotArgs = varargin(2:end);
        end
        plotMode = parseRobustnessPlotMode(plotArgs);
        robustnessFig = plotAggregateRobustnessIndex(robustnessCubeIn,targetMask, ...
            robPassRateIn,method,criteria,slice,ct,cst,plotMode);
        robustnessCube = robustnessFig;
        robPassRate = robPassRateIn;
end
end

function tf = isRobustnessPlotMode(plotMode)
tf = (ischar(plotMode) || (isstring(plotMode) && isscalar(plotMode))) && ...
    any(strcmpi(char(plotMode),{'continuous','binary'}));
end

function tf = isPlotModeName(value)
tf = (ischar(value) || (isstring(value) && isscalar(value))) && ...
    strcmpi(char(value),'plotMode');
end

function plotMode = parseRobustnessPlotMode(args)
plotMode = 'continuous';
if isempty(args)
    return;
end

if isscalar(args) && ~isempty(args{1}) && isRobustnessPlotMode(args{1})
    plotMode = normalizeRobustnessPlotMode(args{1});
    return;
end

for i = 1:2:numel(args)
    if (ischar(args{i}) || (isstring(args{i}) && isscalar(args{i}))) && ...
            strcmpi(char(args{i}),'plotMode') && i < numel(args)
        plotMode = normalizeRobustnessPlotMode(args{i+1});
        return;
    end
end
end

function options = parseRobustnessOptions(args)
options.plotMode = parseRobustnessPlotMode(args);
options.robustnessTargetMode = 'all';
options.robustnessTargets = [];

for i = 1:2:numel(args)
    if ~(ischar(args{i}) || (isstring(args{i}) && isscalar(args{i})))
        continue;
    end

    name = lower(char(args{i}));
    if any(strcmp(name,{'robustnesstargetmode','robustnesstargets'})) && i == numel(args)
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Missing value for robustness option ''%s''.',char(args{i}));
    end

    switch name
        case 'robustnesstargetmode'
            options.robustnessTargetMode = args{i+1};
        case 'robustnesstargets'
            options.robustnessTargets = args{i+1};
    end
end
end

function plotMode = normalizeRobustnessPlotMode(plotMode)
plotMode = lower(char(plotMode));
switch plotMode
    case 'continuous'
        plotMode = 'continuous';
    case 'binary'
        plotMode = 'binary';
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Robustness plot mode must be ''continuous'' or ''binary''.');
end
end

function validateRobustnessInput(meanCube,stdCube,refDose,criteria)
validInput = isequal(size(meanCube),size(stdCube)) && ...
    isnumeric(refDose) && isscalar(refDose) && isfinite(refDose) && refDose > 0 && ...
    isnumeric(criteria) && numel(criteria) == 2 && all(isfinite(criteria)) && all(criteria > 0);
if ~validInput
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Invalid robustness index input.');
end
end

function method = normalizeRobustnessMethod(method)
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

validMethod = ischar(method) && any(strcmp(method,{'index1','index2'}));
if ~validMethod
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Robustness index method must be ''index1'' or ''index2''.');
end
end

function [targetMask,doseMask,refScen] = getRobustnessMasks(meanCube,stdCube,ct,cst,robustnessTargetMode,robustnessTargets)
refScen = getRefScen(ct);
targetMask = NaN(size(meanCube));
doseMask = meanCube > 0 | stdCube > 0;

if isempty(cst)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Robustness index calculation requires a cst input to identify TARGET voxels.');
end

targetRows = matRad_resolveStructureSelection(cst,robustnessTargetMode,robustnessTargets,'TARGET');
for i = targetRows
    hasTarget = isequal(cst{i,3},'TARGET') && numel(cst{i,4}) >= refScen && ...
        ~isempty(cst{i,4}{refScen});
    if hasTarget
        targetMask(cst{i,4}{refScen}) = 1;
    end
end
end

function robPassRate = calcPassRate(passMask,targetMask)
targetIx = ~isnan(targetMask);
numOfTargetVoxels = sum(targetIx(:));
if numOfTargetVoxels > 0
    numOfPassRobustness = sum(passMask(:) & targetIx(:));
    robPassRate = 100 * numOfPassRobustness / numOfTargetVoxels;
else
    robPassRate = NaN;
end
end

function robustnessFig = plotAggregateRobustnessIndex(robustnessCube,targetMask,robPassRate,method,criteria,slice,ct,cst,plotMode)
robustnessFig = [];
if isempty(slice)
    return;
end

meanDoseThreshold = criteria(1);
stdThreshold = criteria(2);
targetMaskCube = NaN(size(robustnessCube));
targetMaskCube(targetMask) = 1;
refScen = getRefScen(ct);

switch method
    case 'index1'
        maxRob = 5.01;
        doseWindow = [0 maxRob];
        mMap1 = round(1/maxRob*256);
        mMap2 = 256 - mMap1;
        colormap1 = [linspace(0.40,1,mMap1)',linspace(1,1,mMap1)',linspace(0.40,1,mMap1)'];
        colormap2 = matRad_getColormap('gammaIndex',2*mMap2);
        myColormap = [colormap1; colormap2(mMap2+1:end-1,:)];

        robustnessFig = figure;
        robustnessFig.Position(3:4) = [400 400];
        if strcmp(plotMode,'binary')
            plotRobustnessSlice(gca,ct,cst,refScen,(robustnessCube < 1).*targetMaskCube, ...
                3,slice,[],[0 2.01],[],'LineWidth',1.5);
        else
            plotRobustnessSlice(gca,ct,cst,refScen,robustnessCube.*targetMaskCube, ...
                3,slice,myColormap,doseWindow,'Delta Index','LineWidth',1.5);
        end
        subtitle(formatRobustnessSubtitleText(robPassRate,'target voxels', ...
            meanDoseThreshold,stdThreshold));
        matRad_addFigureTextbox(robustnessFig,formatRobustnessIndexText(robPassRate), ...
            'tag','matRadMetricTextbox');
    case 'index2'
        robustnessFig = figure;
        robustnessFig.Position(3:4) = [400 400];
        plotRobustnessSlice(gca,ct,cst,refScen,robustnessCube.*targetMaskCube, ...
            3,slice,[],[0 2.01],'Delta Index','LineWidth',1.5);
        subtitle(formatRobustnessSubtitleText(robPassRate,'target voxels', ...
            meanDoseThreshold,stdThreshold));
        matRad_addFigureTextbox(robustnessFig,formatRobustnessIndexText(robPassRate), ...
            'tag','matRadMetricTextbox');
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Robustness index method must be ''index1'' or ''index2''.');
end
end

function robustnessFig = plotRobustnessIndex1(robustnessCube,targetMask,robPassRate,meanDoseThreshold,stdThreshold,slice,ct,cst,refScen,plotMode)
robustnessFig = [];
if isempty(slice)
    return;
end

plane = 3;
maxRob = 5.01;
doseWindow = [0 maxRob];

mMap1 = round(1/maxRob*256);
mMap2 = 256 - mMap1;

colormap1 = [linspace(0.40,1,mMap1)',linspace(1,1,mMap1)',linspace(0.40,1,mMap1)'];
colormap2 = matRad_getColormap('gammaIndex',2*mMap2);
myColormap = [colormap1; colormap2(mMap2+1:end-1,:)];

robustnessFig = figure;
robustnessFig.Position(3:4) = [400 400];

numSlices = getNumSlices(robustnessCube,plane);
if strcmp(plotMode,'binary')
    plotCube = (robustnessCube < 1).*targetMask;
    plotColorMap = [];
    plotWindow = [0 2.01];
    colorBarLabel = [];
else
    plotCube = robustnessCube.*targetMask;
    plotColorMap = myColormap;
    plotWindow = doseWindow;
    colorBarLabel = 'Delta Index';
end

ax1 = gca;
plotRobustnessSlice(ax1,ct,cst,refScen,plotCube, ...
    plane,slice,plotColorMap,plotWindow,colorBarLabel,'LineWidth',1.5);
subtitle(formatRobustnessSubtitleText(robPassRate,'points', ...
    meanDoseThreshold,stdThreshold));
matRad_addFigureTextbox(robustnessFig,formatRobustnessIndexText(robPassRate), ...
    'tag','matRadMetricTextbox');

sliderStep = getSliderStep(numSlices);
b = uicontrol('Parent',robustnessFig,'Style','slider','Position',[70,5,280,23], ...
    'value',slice,'min',1,'max',numSlices,'SliderStep',sliderStep);
b.Callback = @(es,ed) plotRobustnessSlice(ax1,ct,cst,refScen, ...
    plotCube,plane,round(es.Value),plotColorMap,plotWindow,colorBarLabel, ...
    'LineWidth',1.5);
end

function robustnessFig = plotRobustnessIndex2(robustnessCube,targetMask,doseMask,robPassRate,meanDoseThreshold,stdThreshold,slice,ct,cst,refScen)
robustnessFig = [];
if isempty(slice)
    return;
end

plane = 3;
doseWindow = [0 2.01];

robustnessFig = figure;
robustnessFig.Position(3:4) = [400 400];

numSlices = getNumSlices(robustnessCube,plane);
plotCube = robustnessCube.*targetMask;
if all(isnan(plotCube(:)))
    plotCube = robustnessCube.*doseMask;
end
plotRobustnessSlice(gca,ct,cst,refScen,plotCube,plane,slice,[],doseWindow, ...
    'Delta Index','LineWidth',1.5);
subtitle(formatRobustnessSubtitleText(robPassRate,'points', ...
    meanDoseThreshold,stdThreshold));
matRad_addFigureTextbox(robustnessFig,formatRobustnessIndexText(robPassRate), ...
    'tag','matRadMetricTextbox');

ax1 = gca;
sliderStep = getSliderStep(numSlices);
b = uicontrol('Parent',robustnessFig,'Style','slider','Position',[70,5,280,23], ...
    'value',slice,'min',1,'max',numSlices,'SliderStep',sliderStep);
b.Callback = @(es,ed) plotRobustnessSlice(ax1,ct,cst,refScen,plotCube, ...
    plane,round(es.Value),[],doseWindow,'Delta Index','LineWidth',1.5);
end

function refScen = getRefScen(ct)
if isstruct(ct) && isfield(ct,'refScen')
    refScen = ct.refScen;
else
    refScen = 1;
end
end

function plotRobustnessSlice(axesHandle,ct,cst,refScen,plotCube,plane,slice,doseColorMap,doseWindow,colorBarLabel,varargin)
[plotCube,doseColorMap] = matRad_prepareLimitedIndexPlot(plotCube,doseWindow,doseColorMap);
if canUseSliceWrapper(ct)
    matRad_plotSliceWrapper(axesHandle,ct,cst,refScen,plotCube,plane,slice,[],[], ...
        colorcube,doseColorMap,doseWindow,[],[],colorBarLabel,[],varargin{:});
else
    plotCubeSlice(axesHandle,plotCube,plane,slice,doseColorMap,doseWindow,colorBarLabel);
end
end

function plotCubeSlice(axesHandle,plotCube,plane,slice,doseColorMap,doseWindow,colorBarLabel)
if isempty(doseColorMap)
    doseColorMap = jet(64);
end

cla(axesHandle);
imagesc(axesHandle,getCubeSlice(plotCube,plane,slice));
if ~isempty(doseWindow)
    set(axesHandle,'CLim',doseWindow);
end
colormap(axesHandle,doseColorMap);
hCMap = colorbar(axesHandle);
if ~isempty(colorBarLabel)
    set(get(hCMap,'YLabel'),'String',colorBarLabel,'FontSize',14);
end
axis(axesHandle,'image');
set(axesHandle,'xtick',[],'ytick',[]);
end

function cubeSlice = getCubeSlice(cube,plane,slice)
switch plane
    case 1
        cubeSlice = squeeze(cube(slice,:,:));
    case 2
        cubeSlice = squeeze(cube(:,slice,:));
    case 3
        cubeSlice = cube(:,:,slice);
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Invalid plane ''%d'' selected for visualization!',plane);
end
end

function numSlices = getNumSlices(cube,plane)
numSlices = size(cube,plane);
end

function tf = canUseSliceWrapper(ct)
tf = isstruct(ct) && isfield(ct,'cubeHU') && iscell(ct.cubeHU) && ...
    isfield(ct,'resolution') && (isfield(ct,'cubeDim') || isfield(ct,'dimensions'));
end

function sliderStep = getSliderStep(numSlices)
if numSlices > 1
    sliderStep = [1/(numSlices-1) 1/(numSlices-1)];
else
    sliderStep = [1 1];
end
end

function metricText = formatRobustnessIndexText(robPassRate)
if isempty(robPassRate) || ~isnumeric(robPassRate) || ~isscalar(robPassRate) || ~isfinite(robPassRate)
    metricText = 'RI = n/a';
else
    metricText = sprintf('RI = %.4f',robPassRate/100);
end
end

function subtitleText = formatRobustnessSubtitleText(robPassRate,voxelScope,meanDoseThreshold,stdThreshold)
subtitleText = {[num2str(robPassRate,5) '% of ' voxelScope ...
    ' pass robustness criterion (' num2str(meanDoseThreshold) '% / ' ...
    num2str(stdThreshold) '%)']; ' '};
end
