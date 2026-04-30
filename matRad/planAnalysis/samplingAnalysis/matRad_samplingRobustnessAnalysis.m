function [robustnessAnalysis,robustnessFig1,robustnessFig2] = matRad_samplingRobustnessAnalysis(meanCube,stdCube,criteria,ct,cst,pln,slice,varargin)
% matRad_samplingRobustnessAnalysis calculates sampling robustness analysis
%
% call
%   robustnessAnalysis = matRad_samplingRobustnessAnalysis(meanCube,stdCube,criteria,ct,cst,pln)
%   [robustnessAnalysis,robustnessFig1,robustnessFig2] = matRad_samplingRobustnessAnalysis(meanCube,stdCube,criteria,ct,cst,pln,slice)
%
% input
%   meanCube:           expected dose cube
%   stdCube:            dose standard deviation cube
%   criteria:           1x2 vector [mean dose deviation %, std dose %]
%   ct:                 matRad ct struct
%   cst:                matRad cst cell array
%   pln:                matRad pln struct
%   slice:              (optional) CT slice used to create figures
%   varargin:           optional Name/Value pairs:
%                       - 'robustnessTargetMode': 'all', 'include', or
%                         'exclude'
%                       - 'robustnessTargets': target names or cst row
%                         indices for include/exclude mode
%                       - 'sampleMask': logical CT mask marking sampled
%                         voxels. Targets not fully covered are skipped.
%
% output
%   robustnessAnalysis: robustness analysis struct with target-wise and
%                       aggregate robustness indices
%   robustnessFig1:     aggregate robustness figure for index1
%   robustnessFig2:     aggregate robustness figure for index2
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

matRad_cfg = MatRad_Config.instance();
robustnessFig1 = [];
robustnessFig2 = [];

p = inputParser;
p.CaseSensitive = false;
p.addParameter('robustnessTargetMode','all',@(m) ischar(m) || (isstring(m) && isscalar(m)));
p.addParameter('robustnessTargets',[],@(s) isempty(s) || isnumeric(s) || ischar(s) || isstring(s) || iscell(s));
p.addParameter('sampleMask',[],@(m) isempty(m) || (islogical(m) && isequal(size(m),size(meanCube))));
parse(p,varargin{:});

sampleMask = p.Results.sampleMask;
if isempty(sampleMask)
    sampleMask = isfinite(meanCube) & isfinite(stdCube);
end

robustnessAnalysis.meanCubeName = 'doseStat.meanCubeW';
robustnessAnalysis.meanCube = meanCube;
robustnessAnalysis.stdCubeName = 'doseStat.stdCubeW';
robustnessAnalysis.stdCube = stdCube;
robustnessAnalysis.sampleMask = sampleMask;
robustnessAnalysis.sampleCoverageFraction = nnz(sampleMask) / numel(sampleMask);
robustnessAnalysis.criteria = criteria;
robustnessAnalysis.meanDoseThreshold = criteria(1);
robustnessAnalysis.stdThreshold = criteria(2);
robustnessAnalysis.index1.method = 'index1';
robustnessAnalysis.index1.criteria = criteria;
robustnessAnalysis.index2.method = 'index2';
robustnessAnalysis.index2.criteria = criteria;

[selectedTargetRows,structureSelection] = matRad_resolveStructureSelection( ...
    cst,p.Results.robustnessTargetMode,p.Results.robustnessTargets,'TARGET');
targetDoseInfo = matRad_getTargetReferenceDoses(cst,pln);
targetDoseInfo = targetDoseInfo(ismember([targetDoseInfo.cstIndex],selectedTargetRows));
robustnessAnalysis.structureSelection = structureSelection;
robustnessAnalysis.refDose = [targetDoseInfo.refDose];
robustnessAnalysis.targets = targetDoseInfo;

if ~isempty(targetDoseInfo) && any(isfinite([targetDoseInfo.refDose]) & [targetDoseInfo.refDose] > 0)
    matRad_cfg.dispInfo(['matRad: Performing robustness index analysis with parameters ', ...
        num2str(criteria), ' [%% %%] \n']);

    aggregateCube1 = NaN(ct.cubeDim);
    aggregateCube2 = NaN(ct.cubeDim);
    aggregateTargetMask = false(ct.cubeDim);

    for robTarget = 1:numel(targetDoseInfo)
        cstIdx = targetDoseInfo(robTarget).cstIndex;
        refDose = targetDoseInfo(robTarget).refDose;
        targetVoxels = getTargetVoxels(cst,cstIdx,ct);

        robustnessAnalysis.targets(robTarget).numVoxels = numel(targetVoxels);
        robustnessAnalysis.targets(robTarget).numUnsampledVoxels = 0;
        robustnessAnalysis.targets(robTarget).isEvaluable = true;
        robustnessAnalysis.targets(robTarget).samplingStatus = 'evaluable';
        robustnessAnalysis.targets(robTarget).index1.method = 'index1';
        robustnessAnalysis.targets(robTarget).index1.criteria = criteria;
        robustnessAnalysis.targets(robTarget).index2.method = 'index2';
        robustnessAnalysis.targets(robTarget).index2.criteria = criteria;

        if isempty(targetVoxels) || ~isfinite(refDose) || refDose <= 0
            robustnessAnalysis.targets(robTarget) = markTargetNotEvaluable( ...
                robustnessAnalysis.targets(robTarget),'invalidReferenceDoseOrEmptyTarget');
            continue;
        end

        unsampledTargetVoxels = targetVoxels(~sampleMask(targetVoxels));
        if ~isempty(unsampledTargetVoxels)
            robustnessAnalysis.targets(robTarget).numUnsampledVoxels = numel(unsampledTargetVoxels);
            robustnessAnalysis.targets(robTarget) = markTargetNotEvaluable( ...
                robustnessAnalysis.targets(robTarget),'partialSamplingCoverage');
            matRad_cfg.dispWarning(['Skipping robustness target ''%s'' because ', ...
                '%d of %d target voxels were not sampled.\n'], ...
                robustnessAnalysis.targets(robTarget).name,numel(unsampledTargetVoxels),numel(targetVoxels));
            continue;
        end

        targetCst = cst(cstIdx,:);
        [targetRobustnessCube1,targetRobPassRate1] = matRad_robustnessIndex( ...
            meanCube,stdCube,refDose,criteria,'index1',ct,targetCst,[]);
        [targetRobustnessCube2,targetRobPassRate2] = matRad_robustnessIndex( ...
            meanCube,stdCube,refDose,criteria,'index2',ct,targetCst,[]);

        robustnessAnalysis.targets(robTarget).index1.robPassRate = targetRobPassRate1;
        robustnessAnalysis.targets(robTarget).index1.robustnessIndex = targetRobPassRate1/100;
        robustnessAnalysis.targets(robTarget).index1.numPassVoxels = ...
            sum(targetRobustnessCube1(targetVoxels) < 1);
        robustnessAnalysis.targets(robTarget).index2.robPassRate = targetRobPassRate2;
        robustnessAnalysis.targets(robTarget).index2.robustnessIndex = targetRobPassRate2/100;
        robustnessAnalysis.targets(robTarget).index2.numPassVoxels = ...
            sum(targetRobustnessCube2(targetVoxels) == 1);

        aggregateCube1 = updateAggregateRobustnessCube(aggregateCube1, ...
            targetRobustnessCube1,targetVoxels,'index1');
        aggregateCube2 = updateAggregateRobustnessCube(aggregateCube2, ...
            targetRobustnessCube2,targetVoxels,'index2');
        aggregateTargetMask(targetVoxels) = true;
    end

    robustnessAnalysis.evaluableTargetMask = aggregateTargetMask;
    if any(aggregateTargetMask(:))
        robustnessAnalysis.index1.robustnessCube = aggregateCube1;
        robustnessAnalysis.index1.robPassRate = calcAggregatePassRate( ...
            aggregateCube1 < 1,aggregateTargetMask);
        robustnessAnalysis.index1.robustnessIndex = robustnessAnalysis.index1.robPassRate/100;

        robustnessAnalysis.index2.robustnessCube = aggregateCube2;
        robustnessAnalysis.index2.robPassRate = calcAggregatePassRate( ...
            aggregateCube2 == 1,aggregateTargetMask);
        robustnessAnalysis.index2.robustnessIndex = robustnessAnalysis.index2.robPassRate/100;

        if nargout > 1
            robustnessFig1 = matRad_robustnessIndex('plotAggregate',aggregateCube1, ...
                aggregateTargetMask,robustnessAnalysis.index1.robPassRate, ...
                'index1',criteria,slice,ct,cst);
        end
        if nargout > 2
            robustnessFig2 = matRad_robustnessIndex('plotAggregate',aggregateCube2, ...
                aggregateTargetMask,robustnessAnalysis.index2.robPassRate, ...
                'index2',criteria,slice,ct,cst);
        end
    else
        matRad_cfg.dispWarning('No selected target is fully covered by sampled voxels. Skipping aggregate robustness index calculation.\n');
        robustnessAnalysis.index1.robustnessCube = NaN(ct.cubeDim);
        robustnessAnalysis.index1.robPassRate = [];
        robustnessAnalysis.index1.robustnessIndex = [];
        robustnessAnalysis.index2.robustnessCube = NaN(ct.cubeDim);
        robustnessAnalysis.index2.robPassRate = [];
        robustnessAnalysis.index2.robustnessIndex = [];
    end
else
    matRad_cfg.dispWarning('Could not determine target reference dose for any target. Skipping robustness index calculation.\n');
    robustnessAnalysis.index1.robustnessCube = [];
    robustnessAnalysis.index1.robPassRate = [];
    robustnessAnalysis.index1.robustnessIndex = [];
    robustnessAnalysis.index2.robustnessCube = [];
    robustnessAnalysis.index2.robPassRate = [];
    robustnessAnalysis.index2.robustnessIndex = [];
    robustnessAnalysis.evaluableTargetMask = false(ct.cubeDim);
end

end

function targetVoxels = getTargetVoxels(cst,cstIdx,ct)
refScen = getRefScen(ct);
targetVoxels = [];
if numel(cst{cstIdx,4}) >= refScen && ~isempty(cst{cstIdx,4}{refScen})
    targetVoxels = unique(cst{cstIdx,4}{refScen}(:));
end
end

function refScen = getRefScen(ct)
if isfield(ct,'refScen')
    refScen = ct.refScen;
else
    refScen = 1;
end
end

function targetInfo = markTargetNotEvaluable(targetInfo,status)
targetInfo.isEvaluable = false;
targetInfo.samplingStatus = status;
targetInfo.index1.robPassRate = [];
targetInfo.index1.robustnessIndex = [];
targetInfo.index1.numPassVoxels = [];
targetInfo.index2.robPassRate = [];
targetInfo.index2.robustnessIndex = [];
targetInfo.index2.numPassVoxels = [];
end

function aggregateCube = updateAggregateRobustnessCube(aggregateCube,targetCube,targetVoxels,method)
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

function robPassRate = calcAggregatePassRate(passMask,targetMask)
numOfTargetVoxels = sum(targetMask(:));
if numOfTargetVoxels > 0
    robPassRate = 100 * sum(passMask(targetMask)) / numOfTargetVoxels;
else
    robPassRate = NaN;
end
end
