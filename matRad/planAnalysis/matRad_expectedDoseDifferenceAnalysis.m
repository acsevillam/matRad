function expectedDoseDifferenceAnalysis = matRad_expectedDoseDifferenceAnalysis(sampleDoseMatrix, referenceCube, varargin)
% matRad_expectedDoseDifferenceAnalysis calculates expected dose differences
% relative to a reference dose cube
%
% call
%   expectedDoseDifferenceAnalysis = matRad_expectedDoseDifferenceAnalysis(sampleDoseMatrix,referenceCube,...)
%
% input
%   sampleDoseMatrix:   sampled dose matrix with one row per sampled voxel
%                       and one column per sampled scenario
%   referenceCube:      reference dose cube on the same CT grid
%   varargin:           optional Name/Value pairs:
%                       - 'sampledVoxelIndices': linear voxel indices for
%                         the rows of sampleDoseMatrix
%                       - 'sampleMask': logical mask of evaluable voxels
%                       - 'scenWeights': scenario weights
%                       - 'tolerance': dose tolerance around reference
%                       - 'quantity': dose quantity name
%                       - 'evaluationModeBase': internal dose unit mode
%                       - 'referenceName': descriptive name of referenceCube
%
% output
%   expectedDoseDifferenceAnalysis: analysis struct containing status, metadata,
%                        expected dose-difference cubes, auxiliary
%                        over/under/near reference probabilities, and summary values
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

p = inputParser;
p.CaseSensitive = false;
p.addParameter('sampledVoxelIndices', [], @(ix) isempty(ix) || isnumeric(ix));
p.addParameter('sampleMask', [], @(mask) isempty(mask) || islogical(mask));
p.addParameter('scenWeights', [], @(weights) isempty(weights) || isnumeric(weights));
p.addParameter('tolerance', 0, @(value) isnumeric(value) && isscalar(value) && ...
               isfinite(value) && value >= 0);
p.addParameter('quantity', '', @(value) ischar(value) || (isstring(value) && isscalar(value)));
p.addParameter('evaluationModeBase', 'perFraction', @(value) ischar(value) || ...
               (isstring(value) && isscalar(value)));
p.addParameter('referenceName', 'referenceCube', @(value) ischar(value) || ...
               (isstring(value) && isscalar(value)));
parse(p, varargin{:});

expectedDoseDifferenceAnalysis = matRad_initializeExpectedDoseDifferenceAnalysis(referenceCube, p.Results);
[sampledVoxelIndices, sampleMask, isValid, reason] = ...
    matRad_resolveExpectedDoseDifferenceSamplingMask(sampleDoseMatrix, referenceCube, p.Results);
if ~isValid
    expectedDoseDifferenceAnalysis = matRad_markExpectedDoseDifferenceSkipped(expectedDoseDifferenceAnalysis, reason);
    return
end

[weights, isValid, reason] = matRad_resolveExpectedDoseDifferenceWeights(sampleDoseMatrix, p.Results.scenWeights);
if ~isValid
    expectedDoseDifferenceAnalysis = matRad_markExpectedDoseDifferenceSkipped(expectedDoseDifferenceAnalysis, reason);
    return
end

referenceValues = referenceCube(sampledVoxelIndices);
evaluableRows = sampleMask(sampledVoxelIndices) & isfinite(referenceValues) & ...
    all(isfinite(sampleDoseMatrix), 2);
evaluableIndices = sampledVoxelIndices(evaluableRows);

overReferenceProbabilityCube = NaN(size(referenceCube));
underReferenceProbabilityCube = NaN(size(referenceCube));
nearReferenceProbabilityCube = NaN(size(referenceCube));
signedReferenceProbabilityCube = NaN(size(referenceCube));
overReferenceExpectedDoseDifferenceCube = NaN(size(referenceCube));
underReferenceExpectedDoseDifferenceCube = NaN(size(referenceCube));
signedExpectedDoseDifferenceCube = NaN(size(referenceCube));

if any(evaluableRows)
    doseRows = double(sampleDoseMatrix(evaluableRows, :));
    referenceRows = double(referenceValues(evaluableRows));
    tolerance = p.Results.tolerance;
    doseDifferenceRows = doseRows - referenceRows;

    overValues = (doseRows > referenceRows + tolerance) * weights;
    underValues = (doseRows < referenceRows - tolerance) * weights;
    nearValues = (abs(doseRows - referenceRows) <= tolerance) * weights;
    overExpectedDoseDifferenceValues = max(doseDifferenceRows, 0) * weights;
    underExpectedDoseDifferenceValues = max(-doseDifferenceRows, 0) * weights;
    signedExpectedDoseDifferenceValues = doseDifferenceRows * weights;

    overReferenceProbabilityCube(evaluableIndices) = overValues;
    underReferenceProbabilityCube(evaluableIndices) = underValues;
    nearReferenceProbabilityCube(evaluableIndices) = nearValues;
    signedReferenceProbabilityCube(evaluableIndices) = overValues - underValues;
    overReferenceExpectedDoseDifferenceCube(evaluableIndices) = overExpectedDoseDifferenceValues;
    underReferenceExpectedDoseDifferenceCube(evaluableIndices) = underExpectedDoseDifferenceValues;
    signedExpectedDoseDifferenceCube(evaluableIndices) = signedExpectedDoseDifferenceValues;
end

expectedDoseDifferenceAnalysis.status = matRad_resolveExpectedDoseDifferenceStatus(sampleMask);
expectedDoseDifferenceAnalysis.reason = '';
expectedDoseDifferenceAnalysis.sampleMask = sampleMask;
expectedDoseDifferenceAnalysis.sampleCoverageFraction = nnz(sampleMask) / numel(sampleMask);
expectedDoseDifferenceAnalysis.evaluableMask = false(size(referenceCube));
expectedDoseDifferenceAnalysis.evaluableMask(evaluableIndices) = true;
expectedDoseDifferenceAnalysis.evaluableCoverageFraction = nnz(expectedDoseDifferenceAnalysis.evaluableMask) / ...
    numel(expectedDoseDifferenceAnalysis.evaluableMask);
expectedDoseDifferenceAnalysis.weights = weights(:);
expectedDoseDifferenceAnalysis.referenceCube = referenceCube;
expectedDoseDifferenceAnalysis.overReferenceProbabilityCube = overReferenceProbabilityCube;
expectedDoseDifferenceAnalysis.underReferenceProbabilityCube = underReferenceProbabilityCube;
expectedDoseDifferenceAnalysis.nearReferenceProbabilityCube = nearReferenceProbabilityCube;
expectedDoseDifferenceAnalysis.signedReferenceProbabilityCube = signedReferenceProbabilityCube;
expectedDoseDifferenceAnalysis.overReferenceExpectedDoseDifferenceCube = ...
    overReferenceExpectedDoseDifferenceCube;
expectedDoseDifferenceAnalysis.underReferenceExpectedDoseDifferenceCube = ...
    underReferenceExpectedDoseDifferenceCube;
expectedDoseDifferenceAnalysis.signedExpectedDoseDifferenceCube = signedExpectedDoseDifferenceCube;
expectedDoseDifferenceAnalysis.doseWindow = ...
    matRad_getExpectedDoseDifferenceDoseWindow(signedExpectedDoseDifferenceCube);
expectedDoseDifferenceAnalysis.summary = ...
    matRad_summarizeExpectedDoseDifference(overReferenceProbabilityCube, ...
                                           underReferenceProbabilityCube, ...
                                           nearReferenceProbabilityCube, ...
                                           signedReferenceProbabilityCube, ...
                                           overReferenceExpectedDoseDifferenceCube, ...
                                           underReferenceExpectedDoseDifferenceCube, ...
                                           signedExpectedDoseDifferenceCube, ...
                                           expectedDoseDifferenceAnalysis.evaluableMask);

end

function expectedDoseDifferenceAnalysis = matRad_initializeExpectedDoseDifferenceAnalysis(referenceCube, parserResults)
expectedDoseDifferenceAnalysis.status = 'pending';
expectedDoseDifferenceAnalysis.reason = '';
expectedDoseDifferenceAnalysis.quantity = char(parserResults.quantity);
expectedDoseDifferenceAnalysis.evaluationModeBase = char(parserResults.evaluationModeBase);
expectedDoseDifferenceAnalysis.referenceName = char(parserResults.referenceName);
expectedDoseDifferenceAnalysis.tolerance = parserResults.tolerance;
expectedDoseDifferenceAnalysis.sampleMask = false(size(referenceCube));
expectedDoseDifferenceAnalysis.sampleCoverageFraction = 0;
expectedDoseDifferenceAnalysis.evaluableMask = false(size(referenceCube));
expectedDoseDifferenceAnalysis.evaluableCoverageFraction = 0;
expectedDoseDifferenceAnalysis.weights = [];
expectedDoseDifferenceAnalysis.referenceCube = referenceCube;
expectedDoseDifferenceAnalysis.overReferenceProbabilityCube = NaN(size(referenceCube));
expectedDoseDifferenceAnalysis.underReferenceProbabilityCube = NaN(size(referenceCube));
expectedDoseDifferenceAnalysis.nearReferenceProbabilityCube = NaN(size(referenceCube));
expectedDoseDifferenceAnalysis.signedReferenceProbabilityCube = NaN(size(referenceCube));
expectedDoseDifferenceAnalysis.overReferenceExpectedDoseDifferenceCube = NaN(size(referenceCube));
expectedDoseDifferenceAnalysis.underReferenceExpectedDoseDifferenceCube = NaN(size(referenceCube));
expectedDoseDifferenceAnalysis.signedExpectedDoseDifferenceCube = NaN(size(referenceCube));
expectedDoseDifferenceAnalysis.doseWindow = [-1 1];
expectedDoseDifferenceAnalysis.summary = matRad_emptyExpectedDoseDifferenceSummary();
end

function expectedDoseDifferenceAnalysis = matRad_markExpectedDoseDifferenceSkipped(expectedDoseDifferenceAnalysis, reason)
expectedDoseDifferenceAnalysis.status = 'skippedInvalidInput';
expectedDoseDifferenceAnalysis.reason = reason;
end

function [sampledVoxelIndices, sampleMask, isValid, reason] = ...
    matRad_resolveExpectedDoseDifferenceSamplingMask(sampleDoseMatrix, referenceCube, parserResults)
sampledVoxelIndices = parserResults.sampledVoxelIndices(:);
numReferenceVoxels = numel(referenceCube);
numSampledRows = size(sampleDoseMatrix, 1);
isValid = true;
reason = '';

if isempty(sampledVoxelIndices)
    if numSampledRows == numReferenceVoxels
        sampledVoxelIndices = (1:numReferenceVoxels)';
    else
        isValid = false;
        reason = 'sampledVoxelIndices are required when sampleDoseMatrix does not cover the full cube.';
        sampleMask = false(size(referenceCube));
        return
    end
end

if numel(sampledVoxelIndices) ~= numSampledRows
    isValid = false;
    reason = 'Number of sampled voxel indices must match sampleDoseMatrix rows.';
    sampleMask = false(size(referenceCube));
    return
end

if any(sampledVoxelIndices < 1) || any(sampledVoxelIndices > numReferenceVoxels) || ...
        any(round(sampledVoxelIndices) ~= sampledVoxelIndices)
    isValid = false;
    reason = 'sampledVoxelIndices must be valid positive integer voxel indices.';
    sampleMask = false(size(referenceCube));
    return
end

sampleMask = parserResults.sampleMask;
if isempty(sampleMask)
    sampleMask = false(size(referenceCube));
    sampleMask(sampledVoxelIndices) = true;
elseif ~isequal(size(sampleMask), size(referenceCube))
    isValid = false;
    reason = 'sampleMask must match referenceCube dimensions.';
    sampleMask = false(size(referenceCube));
end
end

function [weights, isValid, reason] = matRad_resolveExpectedDoseDifferenceWeights(sampleDoseMatrix, inputWeights)
numScenarios = size(sampleDoseMatrix, 2);
isValid = true;
reason = '';

if isempty(inputWeights)
    weights = ones(numScenarios, 1) ./ numScenarios;
else
    weights = double(inputWeights(:));
end

if numel(weights) ~= numScenarios
    isValid = false;
    reason = 'Number of scenario weights must match sampleDoseMatrix columns.';
    return
end

if any(~isfinite(weights)) || any(weights < 0) || sum(weights) <= 0
    isValid = false;
    reason = 'Scenario weights must be finite, non-negative, and sum to a positive value.';
    return
end

weights = weights ./ sum(weights);
end

function status = matRad_resolveExpectedDoseDifferenceStatus(sampleMask)
if nnz(sampleMask) == numel(sampleMask)
    status = 'computedFullCube';
else
    status = 'computedPartialMask';
end
end

function doseWindow = matRad_getExpectedDoseDifferenceDoseWindow(signedDoseCube)
finiteValues = signedDoseCube(isfinite(signedDoseCube));
if isempty(finiteValues)
    doseWindow = [-1 1];
    return
end

windowAbs = max(abs(finiteValues(:)));
if windowAbs == 0
    windowAbs = 1;
end
doseWindow = [-windowAbs windowAbs];
end

function summary = matRad_summarizeExpectedDoseDifference(overCube, underCube, nearCube, ...
                                                          signedCube, overDoseCube, ...
                                                          underDoseCube, signedDoseCube, ...
                                                          evaluableMask)
summary = matRad_emptyExpectedDoseDifferenceSummary();
overValues = overCube(evaluableMask);
underValues = underCube(evaluableMask);
nearValues = nearCube(evaluableMask);
signedValues = signedCube(evaluableMask);
overDoseValues = overDoseCube(evaluableMask);
underDoseValues = underDoseCube(evaluableMask);
signedDoseValues = signedDoseCube(evaluableMask);

overValues = overValues(isfinite(overValues));
underValues = underValues(isfinite(underValues));
nearValues = nearValues(isfinite(nearValues));
signedValues = signedValues(isfinite(signedValues));
overDoseValues = overDoseValues(isfinite(overDoseValues));
underDoseValues = underDoseValues(isfinite(underDoseValues));
signedDoseValues = signedDoseValues(isfinite(signedDoseValues));

if isempty(overValues)
    return
end

summary.meanOverReferenceProbability = mean(overValues);
summary.maxOverReferenceProbability = max(overValues);
summary.meanUnderReferenceProbability = mean(underValues);
summary.maxUnderReferenceProbability = max(underValues);
summary.meanNearReferenceProbability = mean(nearValues);
summary.maxNearReferenceProbability = max(nearValues);
summary.meanSignedReferenceProbability = mean(signedValues);
summary.maxAbsSignedReferenceProbability = max(abs(signedValues));
summary.meanDirectionalProbability = mean(abs(signedValues));
summary.meanOverReferenceExpectedDoseDifference = mean(overDoseValues);
summary.maxOverReferenceExpectedDoseDifference = max(overDoseValues);
summary.meanUnderReferenceExpectedDoseDifference = mean(underDoseValues);
summary.maxUnderReferenceExpectedDoseDifference = max(underDoseValues);
summary.meanSignedExpectedDoseDifference = mean(signedDoseValues);
summary.maxAbsSignedExpectedDoseDifference = max(abs(signedDoseValues));
end

function summary = matRad_emptyExpectedDoseDifferenceSummary()
summary.meanOverReferenceProbability = NaN;
summary.maxOverReferenceProbability = NaN;
summary.meanUnderReferenceProbability = NaN;
summary.maxUnderReferenceProbability = NaN;
summary.meanNearReferenceProbability = NaN;
summary.maxNearReferenceProbability = NaN;
summary.meanSignedReferenceProbability = NaN;
summary.maxAbsSignedReferenceProbability = NaN;
summary.meanDirectionalProbability = NaN;
summary.meanOverReferenceExpectedDoseDifference = NaN;
summary.maxOverReferenceExpectedDoseDifference = NaN;
summary.meanUnderReferenceExpectedDoseDifference = NaN;
summary.maxUnderReferenceExpectedDoseDifference = NaN;
summary.meanSignedExpectedDoseDifference = NaN;
summary.maxAbsSignedExpectedDoseDifference = NaN;
end
