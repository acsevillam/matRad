function doseStat = matRad_calcSamplingDoseStatistics(ct, pln, mSampDose, scenWeights)
% matRad_calcSamplingDoseStatistics calculates voxel-wise sampling statistics
%
% call
%   doseStat = matRad_calcSamplingDoseStatistics(ct,pln,mSampDose,scenWeights)
%
% input
%   ct:          matRad ct struct
%   pln:         matRad pln struct with subIx
%   mSampDose:   sampled dose matrix with one row per voxel and one column
%                per sampled scenario
%   scenWeights: scenario weights for weighted statistics
%
% output
%   doseStat:    struct with meanCube, stdCube, meanCubeW, stdCubeW,
%                sampleMask, sampledVoxelIndices, and sampleCoverageFraction
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

matRad_cfg = MatRad_Config.instance();

ix = pln.subIx(:);
if size(mSampDose, 1) ~= numel(ix)
    matRad_cfg.dispError('Number of sampled dose rows must match pln.subIx.');
end

if size(mSampDose, 2) ~= numel(scenWeights)
    matRad_cfg.dispError('Number of sampled dose columns must match scenario weights.');
end

doseStat.sampledVoxelIndices = ix;
doseStat.sampleMask = false(ct.cubeDim);
doseStat.sampleMask(ix) = true;
doseStat.sampleCoverageFraction = nnz(doseStat.sampleMask) / numel(doseStat.sampleMask);

doseStat.meanCube = NaN(ct.cubeDim);
doseStat.stdCube = NaN(ct.cubeDim);
doseStat.meanCubeW = NaN(ct.cubeDim);
doseStat.stdCubeW = NaN(ct.cubeDim);

doseStat.meanCube(ix) = mean(mSampDose, 2);
doseStat.stdCube(ix) = std(mSampDose, 1, 2);
doseStat.meanCubeW(ix) = matRad_weightedMean(mSampDose', scenWeights)';
doseStat.stdCubeW(ix) = matRad_weightedStd(mSampDose', scenWeights)';

end

function meanValue = matRad_weightedMean(values, weights)
if exist('weights', 'var') && ~isempty(weights)
    if isvector(values) && isvector(weights)
        meanValue = reshape(weights, 1, []) * reshape(values, [], 1) / sum(weights);
    else
        meanValue = reshape(weights, 1, []) * values ./ sum(weights);
    end
else
    meanValue = mean(values);
end
end

function stdValue = matRad_weightedStd(values, weights)
if exist('weights', 'var') && ~isempty(weights)
    if isvector(values) && isvector(weights)
        values = reshape(values, [], 1);
        weights = reshape(weights, [], 1);
        mu = sum(weights .* values) ./ sum(weights);
        stdValue = sqrt(sum(weights .* (values - mu).^2) ./ sum(weights));
    else
        mu = matRad_weightedMean(values, weights);
        stdValue = sqrt(reshape(weights, 1, []) * ...
                        (bsxfun(@minus, values, mu).^2) ./ sum(weights));
    end
else
    stdValue = std(values, 1);
end
end
