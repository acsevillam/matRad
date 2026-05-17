function [radiusFactor, rank] = matRad_buildDoseIntervalRadiusFactor(centeredScenarioMatrix, ...
                                                                     scenarioWeights, cfg, numBixels)
% matRad_buildDoseIntervalRadiusFactor builds one INTERVAL3 factor
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

if strcmp(cfg.RadiusMode, 'extreme')
    [radiusFactor, rank] = matRad_buildExtremeRadiusFactor(centeredScenarioMatrix, ...
                                                           scenarioWeights, numBixels);
else
    [radiusFactor, rank] = matRad_buildStdRadiusFactor(centeredScenarioMatrix, ...
                                                       scenarioWeights, cfg, numBixels);
end
end

function [radiusFactor, rank] = matRad_buildStdRadiusFactor(centeredScenarioMatrix, ...
                                                            scenarioWeights, cfg, numBixels)
tolerance = 1e-12;

kMax = min(cfg.KMax, numBixels);
if kMax == 0
    [radiusFactor, rank] = matRad_buildZeroRankRadiusFactor(numBixels);
    return
end

positiveWeightIx = scenarioWeights(:) > 0;
scenarioWeights = scenarioWeights(positiveWeightIx);
centeredScenarioMatrix = centeredScenarioMatrix(positiveWeightIx, :);

activeBixels = find(any(centeredScenarioMatrix, 1));
if isempty(activeBixels)
    [radiusFactor, rank] = matRad_buildZeroRankRadiusFactor(numBixels);
    return
end

centeredScenarioMatrix = centeredScenarioMatrix(:, activeBixels);
numScenarios = size(centeredScenarioMatrix, 1);
weightedScenarioMatrix = spdiags(sqrt(scenarioWeights), 0, ...
                                 numScenarios, numScenarios) * centeredScenarioMatrix;

if nnz(weightedScenarioMatrix) == 0 || ...
        norm(weightedScenarioMatrix, 'fro') <= tolerance
    [radiusFactor, rank] = matRad_buildZeroRankRadiusFactor(numBixels);
    return
end

scenarioGram = weightedScenarioMatrix * weightedScenarioMatrix';
scenarioGram = full(0.5 .* (scenarioGram + scenarioGram'));
[leftEigenVectors, eigenValueMatrix] = eig(scenarioGram);
eigenValues = real(diag(eigenValueMatrix));
eigenValues(abs(eigenValues) <= tolerance) = 0;
eigenValues(eigenValues < 0) = 0;
[eigenValues, sortIx] = sort(eigenValues, 'descend');
leftEigenVectors = leftEigenVectors(:, sortIx);

positiveIx = find(eigenValues > tolerance);
if isempty(positiveIx)
    [radiusFactor, rank] = matRad_buildZeroRankRadiusFactor(numBixels);
    return
end

eigenValues = eigenValues(positiveIx);
leftEigenVectors = leftEigenVectors(:, positiveIx);

if strcmp(cfg.KMode, 'dynamic')
    rank = matRad_selectFactorRank(eigenValues, kMax, ...
                                   cfg.RetentionThreshold, tolerance);
else
    rank = min(kMax, numel(eigenValues));
end

if rank == 0
    [radiusFactor, rank] = matRad_buildZeroRankRadiusFactor(numBixels);
    return
end

eigenValues = eigenValues(1:rank);
leftEigenVectors = leftEigenVectors(:, 1:rank);
rightFactorBasis = weightedScenarioMatrix' * leftEigenVectors;
rightFactorBasis = bsxfun(@rdivide, rightFactorBasis, ...
                          sqrt(eigenValues(:))');
activeRadiusFactor = bsxfun(@times, rightFactorBasis, sqrt(eigenValues(:))');

radiusFactor = sparse(double(numBixels), double(rank));
radiusFactor(activeBixels, :) = sparse(activeRadiusFactor);
end

function [radiusFactor, rank] = matRad_buildExtremeRadiusFactor(centeredScenarioMatrix, ...
                                                                scenarioWeights, numBixels)
tolerance = 1e-12;

positiveWeightIx = scenarioWeights(:) > 0;
centeredScenarioMatrix = centeredScenarioMatrix(positiveWeightIx, :);
if isempty(centeredScenarioMatrix)
    [radiusFactor, rank] = matRad_buildZeroRankRadiusFactor(numBixels);
    return
end

delta = max(abs(centeredScenarioMatrix), [], 1);
delta(abs(delta) <= tolerance) = 0;
activeBixels = find(delta);

rank = numel(activeBixels);
if rank == 0
    [radiusFactor, rank] = matRad_buildZeroRankRadiusFactor(numBixels);
    return
end

deltaValues = full(delta(activeBixels(:)));
radiusFactor = sparse(activeBixels(:), (1:rank)', deltaValues(:), ...
                      double(numBixels), double(rank));
end

function rank = matRad_selectFactorRank(eigenValues, kMax, retentionThreshold, tolerance)
maxRank = min(kMax, numel(eigenValues));
if maxRank == 0
    rank = 0;
    return
end

rankEigenValues = eigenValues(1:maxRank);
energyScale = max(abs(rankEigenValues));
if ~isfinite(energyScale) || energyScale <= tolerance
    rank = maxRank;
    return
end

relativeEnergy = (rankEigenValues ./ energyScale).^2;
totalEnergy = sum(relativeEnergy);
if ~isfinite(totalEnergy) || totalEnergy <= 0
    rank = maxRank;
    return
end

thresholdEnergy = retentionThreshold * totalEnergy;
rank = find(cumsum(relativeEnergy) >= thresholdEnergy, 1, 'first');
if isempty(rank)
    rank = maxRank;
end
end

function [radiusFactor, rank] = matRad_buildZeroRankRadiusFactor(numBixels)
rank = 0;
radiusFactor = sparse(numBixels, 0);
end
