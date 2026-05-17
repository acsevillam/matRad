function [sampledDoseColumn, sampleResult] = matRad_calculateSamplingScenario(samplingContext, scenarioId)
% matRad_calculateSamplingScenario calculates one sampled scenario result
%
% call
%   [sampledDoseColumn,sampleResult] = matRad_calculateSamplingScenario(...)
%
% input
%   samplingContext: sampling inputs and reference-grid analysis metadata
%   scenarioId:  canonical scenario id to calculate
%
% output
%   sampledDoseColumn: sampled dose values at subIx as a single column
%   sampleResult:      struct with sample metadata, DVH, and quality
%                      indicators
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

ct = samplingContext.ct;
stf = samplingContext.stf;
cst = samplingContext.cst;
pln = samplingContext.pln;
w = samplingContext.w;
cstEval = samplingContext.cstEval;
subIx = samplingContext.subIx;
dvhPoints = samplingContext.dvhPoints;
refGy = samplingContext.refGy;
refVol = samplingContext.refVol;
doseMapping = samplingContext.doseMapping;
quantity = samplingContext.quantity;

scenario = pln.multScen.getScenario(scenarioId);
plnSamp = pln;
plnSamp.multScen = pln.multScen.extractSingleScenario(scenarioId);

resultSamp = matRad_calcDoseForward(ct, cst, stf, plnSamp, w);
doseCubeSample = resultSamp.(quantity);
sampleCtScenId = scenario.ctScenId;

[doseCubeSample, doseMappingInfo] = matRad_mapSampleDoseToReferenceIfNeeded(ct, doseCubeSample, ...
                                                                            sampleCtScenId, doseMapping);

sampledDoseColumn = single(reshape(doseCubeSample(subIx), [], 1));
sampleResult = matRad_createEmptySamplingResult();
sampleResult.scenarioId = scenarioId;
sampleResult.ctScenId = sampleCtScenId;
sampleResult.refScen = doseMapping.refScen;
sampleResult.scenario = scenario;
sampleResult.doseMapping = doseMappingInfo;
sampleResult.analysisQuantity = quantity;
sampleResult.evaluationModeBase = 'perFraction';
sampleResult.dvh = matRad_calcDVH(cstEval, doseCubeSample, 'cum', dvhPoints);
sampleResult.qi = matRad_calcQualityIndicators(cstEval, pln, doseCubeSample, refGy, refVol);

end

function sampleResult = matRad_createEmptySamplingResult()
sampleResult = struct();
sampleResult.scenarioId = [];
sampleResult.ctScenId = [];
sampleResult.refScen = [];
sampleResult.scenario = [];
sampleResult.doseMapping = [];
sampleResult.analysisQuantity = [];
sampleResult.evaluationModeBase = 'perFraction';
sampleResult.dvh = [];
sampleResult.qi = [];
end

function [doseCubeSample, doseMappingInfo] = matRad_mapSampleDoseToReferenceIfNeeded(ct, doseCubeSample, ctScenId, doseMapping)
if doseMapping.enabled
    [doseCubeSample, doseMappingInfo] = matRad_mapDoseToReferenceScenario(ct, doseCubeSample, ctScenId, ...
                                                                          doseMapping.refScen);
else
    doseMappingInfo = struct();
    doseMappingInfo.sourceCtScenId = ctScenId;
    doseMappingInfo.referenceCtScenId = doseMapping.refScen;
    doseMappingInfo.method = doseMapping.method;
    doseMappingInfo.mapped = false;
end
end
