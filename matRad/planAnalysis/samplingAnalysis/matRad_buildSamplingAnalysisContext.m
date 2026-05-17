function analysisContext = matRad_buildSamplingAnalysisContext(ct, pln, caSampRes, mSampDose, doseStat, quantity, scenWeights)
% matRad_buildSamplingAnalysisContext validates common sampling analysis metadata
%
% call
%   analysisContext = matRad_buildSamplingAnalysisContext(ct,pln,caSampRes,mSampDose,doseStat,quantity,scenWeights)
%
% input
%   ct:          matRad ct struct
%   pln:         matRad pln struct with subIx
%   caSampRes:   sampling result struct array
%   mSampDose:   sampled dose matrix on the reference CT grid
%   doseStat:    sampling dose statistics for mSampDose
%   quantity:    analysis quantity
%   scenWeights: scenario weights matching mSampDose columns
%
% output
%   analysisContext: validated common context for all sampling analyses
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

matRadCfg = MatRad_Config.instance();

if iscell(caSampRes)
    caSampRes = [caSampRes{:}];
end

quantity = matRad_validateSamplingAnalysisQuantity(quantity, matRadCfg);
referenceCtScen = matRad_resolveSamplingAnalysisReferenceCt(ct, matRadCfg);
numSamples = matRad_validateSamplingAnalysisSampleCount(caSampRes, mSampDose, scenWeights, matRadCfg);
[sampledVoxelIndices, sampleMask, sampleCoverageFraction] = ...
    matRad_validateSamplingAnalysisVoxelContext(ct, pln, mSampDose, doseStat, matRadCfg);
[scenarioIds, ctScenIds, doseMapping] = ...
    matRad_validateSamplingAnalysisScenarioContext(caSampRes, quantity, referenceCtScen, matRadCfg);

analysisContext = struct();
analysisContext.analysisGrid = 'referenceCt';
analysisContext.referenceCtScen = referenceCtScen;
analysisContext.quantity = quantity;
analysisContext.evaluationModeBase = 'perFraction';
analysisContext.sampledVoxelIndices = sampledVoxelIndices;
analysisContext.sampleMask = sampleMask;
analysisContext.sampleCoverageFraction = sampleCoverageFraction;
analysisContext.scenarioIds = scenarioIds;
analysisContext.ctScenIds = ctScenIds;
analysisContext.scenWeights = scenWeights(:);
analysisContext.doseMapping = doseMapping;
analysisContext.numSamples = numSamples;

end

function quantity = matRad_validateSamplingAnalysisQuantity(quantity, matRadCfg)
if ~(ischar(quantity) || (isstring(quantity) && isscalar(quantity))) || isempty(char(quantity))
    matRadCfg.dispError('Sampling analysis quantity must be a non-empty scalar text value.');
end
quantity = char(quantity);
end

function referenceCtScen = matRad_resolveSamplingAnalysisReferenceCt(ct, matRadCfg)
if isfield(ct, 'refScen') && ~isempty(ct.refScen)
    referenceCtScen = ct.refScen;
else
    referenceCtScen = 1;
end

if ~(isnumeric(referenceCtScen) && isscalar(referenceCtScen) && isfinite(referenceCtScen) && ...
     round(referenceCtScen) == referenceCtScen && referenceCtScen >= 1)
    matRadCfg.dispError('ct.refScen must be a positive integer scalar when provided.');
end
referenceCtScen = double(referenceCtScen);
end

function numSamples = matRad_validateSamplingAnalysisSampleCount(caSampRes, mSampDose, scenWeights, matRadCfg)
numSamples = size(mSampDose, 2);

if numSamples < 1
    matRadCfg.dispError('Sampling analysis requires at least one sampled scenario.');
end

if numel(caSampRes) ~= numSamples
    matRadCfg.dispError('Number of sampling results must match sampled dose columns.');
end

if numel(scenWeights) ~= numSamples
    matRadCfg.dispError('Number of scenario weights must match sampled dose columns.');
end
end

function [sampledVoxelIndices, sampleMask, sampleCoverageFraction] = ...
    matRad_validateSamplingAnalysisVoxelContext(ct, pln, mSampDose, doseStat, matRadCfg)
if ~isfield(pln, 'subIx') || isempty(pln.subIx)
    matRadCfg.dispError('Sampling analysis requires pln.subIx.');
end

sampledVoxelIndices = pln.subIx(:);
numVoxels = prod(ct.cubeDim);

if size(mSampDose, 1) ~= numel(sampledVoxelIndices)
    matRadCfg.dispError('Number of sampled dose rows must match pln.subIx.');
end

if any(~isfinite(sampledVoxelIndices)) || any(round(sampledVoxelIndices) ~= sampledVoxelIndices) || ...
        any(sampledVoxelIndices < 1) || any(sampledVoxelIndices > numVoxels)
    matRadCfg.dispError('pln.subIx must contain valid positive integer voxel indices.');
end

if numel(unique(sampledVoxelIndices)) ~= numel(sampledVoxelIndices)
    matRadCfg.dispError('pln.subIx must not contain duplicate voxel indices.');
end

if ~isfield(doseStat, 'sampledVoxelIndices') || ...
        ~isequal(doseStat.sampledVoxelIndices(:), sampledVoxelIndices)
    matRadCfg.dispError('doseStat.sampledVoxelIndices must match pln.subIx.');
end

if ~isfield(doseStat, 'sampleMask') || ~islogical(doseStat.sampleMask) || ...
        ~matRad_cubeSizeMatches(doseStat.sampleMask, ct.cubeDim)
    matRadCfg.dispError('doseStat.sampleMask must be a logical cube matching ct.cubeDim.');
end

expectedMask = false(ct.cubeDim);
expectedMask(sampledVoxelIndices) = true;
if ~isequal(doseStat.sampleMask, expectedMask)
    matRadCfg.dispError('doseStat.sampleMask must match pln.subIx on the CT reference grid.');
end

sampleMask = doseStat.sampleMask;
sampleCoverageFraction = nnz(sampleMask) / numel(sampleMask);
if ~isfield(doseStat, 'sampleCoverageFraction') || ...
        abs(doseStat.sampleCoverageFraction - sampleCoverageFraction) > sqrt(eps)
    matRadCfg.dispError('doseStat.sampleCoverageFraction is inconsistent with doseStat.sampleMask.');
end
end

function tf = matRad_cubeSizeMatches(value, cubeDim)
valueSize = size(value);
if numel(valueSize) < numel(cubeDim)
    valueSize = [valueSize ones(1, numel(cubeDim) - numel(valueSize))];
end
tf = isequal(valueSize, cubeDim);
end

function [scenarioIds, ctScenIds, doseMapping] = ...
    matRad_validateSamplingAnalysisScenarioContext(caSampRes, quantity, referenceCtScen, matRadCfg)
numSamples = numel(caSampRes);
scenarioIds = zeros(numSamples, 1);
ctScenIds = zeros(numSamples, 1);

if ~isfield(caSampRes, 'doseMapping') || isempty(caSampRes(1).doseMapping)
    matRadCfg.dispError('Sampling results must contain doseMapping metadata.');
end
doseMapping = repmat(caSampRes(1).doseMapping, numSamples, 1);

for i = 1:numSamples
    sample = caSampRes(i);
    scenarioIds(i) = matRad_getPositiveIntegerField(sample, 'scenarioId', i, matRadCfg);
    ctScenIds(i) = matRad_getPositiveIntegerField(sample, 'ctScenId', i, matRadCfg);
    matRad_validateSamplingResultReference(sample, referenceCtScen, i, matRadCfg);
    matRad_validateSamplingResultQuantity(sample, quantity, i, matRadCfg);
    matRad_validateSamplingResultEvaluationMode(sample, i, matRadCfg);
    doseMapping(i) = matRad_validateSamplingResultDoseMapping(sample, ctScenIds(i), referenceCtScen, ...
                                                              i, matRadCfg);
end
end

function value = matRad_getPositiveIntegerField(sample, fieldName, sampleIx, matRadCfg)
if ~isfield(sample, fieldName) || isempty(sample.(fieldName))
    matRadCfg.dispError('Sampling result %d is missing %s.', sampleIx, fieldName);
end

value = sample.(fieldName);
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && round(value) == value && value >= 1)
    matRadCfg.dispError('Sampling result %d field %s must be a positive integer scalar.', ...
                        sampleIx, fieldName);
end
value = double(value);
end

function matRad_validateSamplingResultReference(sample, referenceCtScen, sampleIx, matRadCfg)
if ~isfield(sample, 'refScen') || isempty(sample.refScen)
    matRadCfg.dispError('Sampling result %d is missing refScen.', sampleIx);
end

if sample.refScen ~= referenceCtScen
    matRadCfg.dispError('Sampling result %d refScen does not match the reference CT scenario.', sampleIx);
end
end

function matRad_validateSamplingResultQuantity(sample, quantity, sampleIx, matRadCfg)
if ~isfield(sample, 'analysisQuantity') || isempty(sample.analysisQuantity)
    matRadCfg.dispError('Sampling result %d is missing analysisQuantity.', sampleIx);
end

if ~strcmp(char(sample.analysisQuantity), quantity)
    matRadCfg.dispError('Sampling result %d analysisQuantity does not match the analysis quantity.', ...
                        sampleIx);
end
end

function matRad_validateSamplingResultEvaluationMode(sample, sampleIx, matRadCfg)
if ~isfield(sample, 'evaluationModeBase') || isempty(sample.evaluationModeBase)
    matRadCfg.dispError('Sampling result %d is missing evaluationModeBase.', sampleIx);
end

if ~strcmp(char(sample.evaluationModeBase), 'perFraction')
    matRadCfg.dispError('Sampling result %d evaluationModeBase must be perFraction.', sampleIx);
end
end

function doseMapping = ...
    matRad_validateSamplingResultDoseMapping(sample, ctScenId, referenceCtScen, sampleIx, matRadCfg)
if ~isfield(sample, 'doseMapping') || isempty(sample.doseMapping)
    matRadCfg.dispError('Sampling result %d is missing doseMapping metadata.', sampleIx);
end

doseMapping = sample.doseMapping;
requiredFields = {'sourceCtScenId', 'referenceCtScenId', 'mapped'};
for fieldIx = 1:numel(requiredFields)
    if ~isfield(doseMapping, requiredFields{fieldIx})
        matRadCfg.dispError('Sampling result %d doseMapping is missing %s.', ...
                            sampleIx, requiredFields{fieldIx});
    end
end

if doseMapping.sourceCtScenId ~= ctScenId
    matRadCfg.dispError('Sampling result %d doseMapping source CT scenario does not match ctScenId.', ...
                        sampleIx);
end

if doseMapping.referenceCtScenId ~= referenceCtScen
    matRadCfg.dispError(['Sampling result %d doseMapping reference CT scenario does not match ', ...
                         'the analysis reference CT scenario.'], sampleIx);
end

if ctScenId ~= referenceCtScen && ~doseMapping.mapped
    matRadCfg.dispError(['Sampling result %d is from a non-reference CT scenario but was not ', ...
                         'mapped to the reference CT.'], sampleIx);
end
end
