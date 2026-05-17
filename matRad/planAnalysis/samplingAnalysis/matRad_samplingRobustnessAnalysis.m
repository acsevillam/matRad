function [robustnessAnalysis, robustnessFig1, robustnessFig2] = ...
    matRad_samplingRobustnessAnalysis(meanCube, stdCube, criteria, ct, cst, pln, slice, varargin)
% matRad_samplingRobustnessAnalysis calculates target-wise sampling robustness
%
% call
%   robustnessAnalysis = matRad_samplingRobustnessAnalysis(meanCube,stdCube,criteria,ct,cst,pln)
%   [robustnessAnalysis,robustnessFig1,robustnessFig2] = matRad_samplingRobustnessAnalysis(...)
%
% input
%   meanCube:           expected dose cube in per-fraction units
%   stdCube:            dose standard deviation cube in per-fraction units
%   criteria:           1x2 vector [mean dose deviation %, std dose %]
%   ct:                 matRad ct struct
%   cst:                matRad cst cell array
%   pln:                matRad plan struct
%   slice:              optional CT slice used for figures
%   varargin:           optional Name/Value pairs:
%                       - 'robustnessTargetMode': 'all', 'include', or
%                         'exclude'
%                       - 'robustnessTargets': target names or cst row
%                         indices
%                       - 'sampleMask': logical CT mask of sampled voxels
%
% output
%   robustnessAnalysis: target-wise and aggregate robustness analysis
%   robustnessFig1:     optional index1 figure
%   robustnessFig2:     optional index2 figure
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

if nargin < 7
    slice = [];
end

matRadCfg = MatRad_Config.instance();
robustnessFig1 = [];
robustnessFig2 = [];

p = inputParser;
p.CaseSensitive = false;
p.addParameter('robustnessTargetMode', 'all', @(mode) ischar(mode) || ...
               (isstring(mode) && isscalar(mode)));
p.addParameter('robustnessTargets', [], @(selection) isempty(selection) || ...
               isnumeric(selection) || ischar(selection) || isstring(selection) || ...
               iscell(selection));
p.addParameter('sampleMask', [], @(mask) isempty(mask) || ...
               (islogical(mask) && isequal(size(mask), size(meanCube))));
parse(p, varargin{:});

sampleMask = p.Results.sampleMask;
if isempty(sampleMask)
    sampleMask = isfinite(meanCube) & isfinite(stdCube);
end

robustnessAnalysis = matRad_initializeRobustnessAnalysis(meanCube, stdCube, criteria, ...
                                                         sampleMask, p.Results);

[targetRows, structureSelection] = matRad_resolveStructureSelection( ...
                                                                    cst, p.Results.robustnessTargetMode, p.Results.robustnessTargets, 'TARGET');
targetDoseInfo = matRad_getTargetReferenceDoses(cst, pln);
targetDoseInfo = targetDoseInfo(ismember([targetDoseInfo.cstIndex], targetRows));

robustnessAnalysis.structureSelection = structureSelection;
robustnessAnalysis.refDose = [targetDoseInfo.refDose];
robustnessAnalysis.targets = matRad_initializeTargetRobustness(targetDoseInfo);

if isempty(targetDoseInfo)
    matRadCfg.dispWarning('No target selected for robustness index analysis.\n');
    robustnessAnalysis = matRad_setEmptyAggregate(robustnessAnalysis, ct);
    return
end

matRadCfg.dispInfo(['matRad: Performing robustness index analysis with parameters ', ...
                    num2str(criteria), ' [%% %%] \n']);

[aggregateCube1, aggregateCube2, aggregateTargetMask] = matRad_initializeAggregate(ct);

for targetIx = 1:numel(targetDoseInfo)
    cstIdx = targetDoseInfo(targetIx).cstIndex;
    refDose = targetDoseInfo(targetIx).refDose;
    targetVoxels = matRad_getTargetVoxels(cst, cstIdx, ct);

    robustnessAnalysis.targets(targetIx).numVoxels = numel(targetVoxels);
    robustnessAnalysis.targets(targetIx).numUnsampledVoxels = 0;
    robustnessAnalysis.targets(targetIx).isEvaluable = true;
    robustnessAnalysis.targets(targetIx).samplingStatus = 'evaluable';

    if isempty(targetVoxels) || ~isfinite(refDose) || refDose <= 0
        robustnessAnalysis.targets(targetIx) = matRad_markTargetNotEvaluable( ...
                                                                             robustnessAnalysis.targets(targetIx), ...
                                                                             'invalidReferenceDoseOrEmptyTarget');
        continue
    end

    unsampledTargetVoxels = targetVoxels(~sampleMask(targetVoxels));
    if ~isempty(unsampledTargetVoxels)
        robustnessAnalysis.targets(targetIx).numUnsampledVoxels = numel(unsampledTargetVoxels);
        robustnessAnalysis.targets(targetIx) = matRad_markTargetNotEvaluable( ...
                                                                             robustnessAnalysis.targets(targetIx), 'partialSamplingCoverage');
        matRadCfg.dispWarning(['Skipping robustness target ''%s'' because ', ...
                               '%d of %d target voxels were not sampled.\n'], ...
                              robustnessAnalysis.targets(targetIx).name, ...
                              numel(unsampledTargetVoxels), numel(targetVoxels));
        continue
    end

    targetCst = cst(cstIdx, :);
    [targetCube1, targetPassRate1] = matRad_robustnessIndex( ...
                                                            meanCube, stdCube, refDose, criteria, 'index1', ct, targetCst, []);
    [targetCube2, targetPassRate2] = matRad_robustnessIndex( ...
                                                            meanCube, stdCube, refDose, criteria, 'index2', ct, targetCst, []);

    robustnessAnalysis.targets(targetIx).index1.robPassRate = targetPassRate1;
    robustnessAnalysis.targets(targetIx).index1.robustnessIndex = targetPassRate1 / 100;
    robustnessAnalysis.targets(targetIx).index1.numPassVoxels = ...
        sum(targetCube1(targetVoxels) < 1);

    robustnessAnalysis.targets(targetIx).index2.robPassRate = targetPassRate2;
    robustnessAnalysis.targets(targetIx).index2.robustnessIndex = targetPassRate2 / 100;
    robustnessAnalysis.targets(targetIx).index2.numPassVoxels = ...
        sum(targetCube2(targetVoxels) == 1);

    aggregateCube1 = matRad_updateAggregateRobustnessCube(aggregateCube1, ...
                                                          targetCube1, targetVoxels, 'index1');
    aggregateCube2 = matRad_updateAggregateRobustnessCube(aggregateCube2, ...
                                                          targetCube2, targetVoxels, 'index2');
    aggregateTargetMask(targetVoxels) = true;
end

robustnessAnalysis.evaluableTargetMask = aggregateTargetMask;
robustnessAnalysis = matRad_setAggregateResults(robustnessAnalysis, aggregateCube1, ...
                                                aggregateCube2, aggregateTargetMask, ct);

if ~isempty(slice) && any(aggregateTargetMask(:))
    [robustnessFig1, robustnessFig2] = matRad_plotAggregateRobustness( ...
                                                                      robustnessAnalysis, aggregateTargetMask, criteria, ct, cst, slice);
end

end

function robustnessAnalysis = matRad_initializeRobustnessAnalysis(meanCube, stdCube, criteria, sampleMask, parserResults)
robustnessAnalysis.meanCubeName = 'doseStat.meanCubeW';
robustnessAnalysis.meanCube = meanCube;
robustnessAnalysis.stdCubeName = 'doseStat.stdCubeW';
robustnessAnalysis.stdCube = stdCube;
robustnessAnalysis.sampleMask = sampleMask;
robustnessAnalysis.sampleCoverageFraction = nnz(sampleMask) / numel(sampleMask);
robustnessAnalysis.criteria = criteria;
robustnessAnalysis.meanDoseThreshold = criteria(1);
robustnessAnalysis.stdThreshold = criteria(2);
robustnessAnalysis.robustnessTargetMode = parserResults.robustnessTargetMode;
robustnessAnalysis.robustnessTargets = parserResults.robustnessTargets;
robustnessAnalysis.index1.method = 'index1';
robustnessAnalysis.index1.criteria = criteria;
robustnessAnalysis.index2.method = 'index2';
robustnessAnalysis.index2.criteria = criteria;
end

function targets = matRad_initializeTargetRobustness(targetDoseInfo)
targets = targetDoseInfo;
for targetIx = 1:numel(targets)
    targets(targetIx).numVoxels = [];
    targets(targetIx).numUnsampledVoxels = [];
    targets(targetIx).isEvaluable = [];
    targets(targetIx).samplingStatus = '';
    targets(targetIx).index1.method = 'index1';
    targets(targetIx).index1.criteria = [];
    targets(targetIx).index1.robPassRate = [];
    targets(targetIx).index1.robustnessIndex = [];
    targets(targetIx).index1.numPassVoxels = [];
    targets(targetIx).index2.method = 'index2';
    targets(targetIx).index2.criteria = [];
    targets(targetIx).index2.robPassRate = [];
    targets(targetIx).index2.robustnessIndex = [];
    targets(targetIx).index2.numPassVoxels = [];
end
end

function [aggregateCube1, aggregateCube2, aggregateTargetMask] = matRad_initializeAggregate(ct)
aggregateCube1 = NaN(ct.cubeDim);
aggregateCube2 = NaN(ct.cubeDim);
aggregateTargetMask = false(ct.cubeDim);
end

function targetVoxels = matRad_getTargetVoxels(cst, cstIdx, ct)
refScen = matRad_getRefScen(ct);
targetVoxels = [];
if numel(cst{cstIdx, 4}) >= refScen && ~isempty(cst{cstIdx, 4}{refScen})
    targetVoxels = unique(cst{cstIdx, 4}{refScen}(:));
end
end

function refScen = matRad_getRefScen(ct)
if isfield(ct, 'refScen')
    refScen = ct.refScen;
else
    refScen = 1;
end
end

function targetInfo = matRad_markTargetNotEvaluable(targetInfo, status)
targetInfo.isEvaluable = false;
targetInfo.samplingStatus = status;
targetInfo.index1.robPassRate = [];
targetInfo.index1.robustnessIndex = [];
targetInfo.index1.numPassVoxels = [];
targetInfo.index2.robPassRate = [];
targetInfo.index2.robustnessIndex = [];
targetInfo.index2.numPassVoxels = [];
end

function aggregateCube = matRad_updateAggregateRobustnessCube(aggregateCube, targetCube, targetVoxels, method)
targetValues = double(targetCube(targetVoxels));
aggregateValues = aggregateCube(targetVoxels);
newVoxels = isnan(aggregateValues);
aggregateValues(newVoxels) = targetValues(newVoxels);

overlapVoxels = ~newVoxels;
if any(overlapVoxels)
    switch method
        case 'index1'
            aggregateValues(overlapVoxels) = max(aggregateValues(overlapVoxels), ...
                                                 targetValues(overlapVoxels));
        case 'index2'
            aggregateValues(overlapVoxels) = min(aggregateValues(overlapVoxels), ...
                                                 targetValues(overlapVoxels));
    end
end

aggregateCube(targetVoxels) = aggregateValues;
end

function robustnessAnalysis = matRad_setAggregateResults(robustnessAnalysis, aggregateCube1, aggregateCube2, aggregateTargetMask, ct)
if any(aggregateTargetMask(:))
    robustnessAnalysis.index1.robustnessCube = aggregateCube1;
    robustnessAnalysis.index1.robPassRate = matRad_calcAggregatePassRate(aggregateCube1 < 1, ...
                                                                         aggregateTargetMask);
    robustnessAnalysis.index1.robustnessIndex = robustnessAnalysis.index1.robPassRate / 100;

    robustnessAnalysis.index2.robustnessCube = aggregateCube2;
    robustnessAnalysis.index2.robPassRate = matRad_calcAggregatePassRate(aggregateCube2 == 1, ...
                                                                         aggregateTargetMask);
    robustnessAnalysis.index2.robustnessIndex = robustnessAnalysis.index2.robPassRate / 100;
else
    matRadCfg = MatRad_Config.instance();
    matRadCfg.dispWarning(['No selected target is fully covered by sampled voxels. ', ...
                           'Skipping aggregate robustness index calculation.\n']);
    robustnessAnalysis = matRad_setEmptyAggregate(robustnessAnalysis, ct);
end
end

function robustnessAnalysis = matRad_setEmptyAggregate(robustnessAnalysis, ct)
robustnessAnalysis.evaluableTargetMask = false(ct.cubeDim);
robustnessAnalysis.index1.robustnessCube = NaN(ct.cubeDim);
robustnessAnalysis.index1.robPassRate = [];
robustnessAnalysis.index1.robustnessIndex = [];
robustnessAnalysis.index2.robustnessCube = NaN(ct.cubeDim);
robustnessAnalysis.index2.robPassRate = [];
robustnessAnalysis.index2.robustnessIndex = [];
end

function robPassRate = matRad_calcAggregatePassRate(passMask, targetMask)
numTargetVoxels = nnz(targetMask);
if numTargetVoxels > 0
    robPassRate = 100 * nnz(passMask(targetMask)) / numTargetVoxels;
else
    robPassRate = NaN;
end
end

function [robustnessFig1, robustnessFig2] = matRad_plotAggregateRobustness(robustnessAnalysis, targetMask, criteria, ct, cst, slice)
robustnessFig1 = matRad_plotAggregateCube(robustnessAnalysis.index1.robustnessCube, ...
                                          targetMask, robustnessAnalysis.index1.robPassRate, ...
                                          criteria, 'index1', ct, cst, slice);
robustnessFig2 = matRad_plotAggregateCube(robustnessAnalysis.index2.robustnessCube, ...
                                          targetMask, robustnessAnalysis.index2.robPassRate, ...
                                          criteria, 'index2', ct, cst, slice);
end

function robustnessFig = matRad_plotAggregateCube(aggregateCube, targetMask, robPassRate, criteria, method, ct, cst, slice)
plotCube = double(aggregateCube);
plotCube(~targetMask) = NaN;

targetCst = cst;
refScen = matRad_getRefScen(ct);
for cstIx = 1:size(targetCst, 1)
    if numel(targetCst{cstIx, 4}) >= refScen
        targetCst{cstIx, 4}{1} = targetCst{cstIx, 4}{refScen};
    end
end

robustnessFig = figure;
set(robustnessFig, 'Color', [1 1 1]);
matRad_plotSlice(ct, 'axesHandle', gca, 'cst', targetCst, 'cubeIdx', 1, ...
                 'dose', plotCube, 'plane', 3, 'slice', slice, ...
                 'contourColorMap', colorcube, ...
                 'doseWindow', matRad_getDoseWindow(method), ...
                 'doseColorMap', matRad_getColorMap(method), ...
                 'colorBarLabel', matRad_getColorbarLabel(method));
title(sprintf('RI = %.4f (%g%% / %g%%)', robPassRate / 100, criteria(1), criteria(2)));
end

function colorMap = matRad_getColorMap(method)
switch method
    case 'index1'
        colorMap = matRad_getColormap('gammaIndex');
    otherwise
        colorMap = [];
end
end

function doseWindow = matRad_getDoseWindow(method)
switch method
    case 'index1'
        doseWindow = [0 5];
    otherwise
        doseWindow = [0 2];
end
end

function label = matRad_getColorbarLabel(method)
switch method
    case 'index1'
        label = 'Delta Index';
    otherwise
        label = [];
end
end
