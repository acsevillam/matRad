function [sampledDoseColumn,sampleResult] = matRad_calculateSamplingScenario(ct,stf,cst,pln,w,cstEval,subIx,dvhPoints,refGy,refVol,samplingCtScenIds,doseMapping,refScen,sampleIx)
% matRad_calculateSamplingScenario calculates one sampled scenario result
%
% call
%   [sampledDoseColumn,sampleResult] = matRad_calculateSamplingScenario(...)
%
% input
%   ct:               matRad ct struct
%   stf:              matRad steering information struct
%   cst:              matRad cst struct
%   pln:              matRad plan meta information struct
%   w:                optional bixel weight vector
%   cstEval:          cst struct in the reference CT scenario
%   subIx:            sampled voxel indices
%   dvhPoints:        common DVH dose grid
%   refGy:            reference dose points for quality indicators
%   refVol:           reference volume points for quality indicators
%   samplingCtScenIds: CT scenario id of each sample
%   doseMapping:      sampling dose mapping configuration
%   refScen:          reference CT scenario id
%   sampleIx:         sample index to calculate
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

plnSamp          = pln;
plnSamp.multScen = pln.multScen.extractSingleScenario(sampleIx);

resultSamp     = matRad_calcDoseForward(ct,cst,stf,plnSamp,w);
doseCubeSample = resultSamp.(pln.bioParam.quantityVis);
sampleCtScenId = samplingCtScenIds(sampleIx);

[doseCubeSample,doseMappingInfo] = mapSampleDoseToReferenceIfNeeded( ...
    ct,doseCubeSample,sampleCtScenId,doseMapping);

sampledDoseColumn = single(reshape(doseCubeSample(subIx),[],1));
sampleResult = createEmptySamplingResult();
sampleResult.bioParam      = pln.bioParam;
sampleResult.ctScenId      = sampleCtScenId;
sampleResult.refScen       = refScen;
sampleResult.doseMapping   = doseMappingInfo;
rangeShift = plnSamp.multScen.getRangeShift(1);
sampleResult.relRangeShift = rangeShift(2);
sampleResult.absRangeShift = rangeShift(1);
sampleResult.isoShift      = plnSamp.multScen.getSetupShift(1);
sampleResult.evaluationModeBase = 'perFraction';
sampleResult.dvh           = matRad_calcDVH(cstEval,doseCubeSample,'cum',dvhPoints);
sampleResult.qi            = matRad_calcQualityIndicators(cstEval,pln,doseCubeSample,refGy,refVol);

end

function sampleResult = createEmptySamplingResult()
sampleResult = struct( ...
    'bioParam',[], ...
    'ctScenId',[], ...
    'refScen',[], ...
    'doseMapping',[], ...
    'relRangeShift',[], ...
    'absRangeShift',[], ...
    'isoShift',[], ...
    'evaluationModeBase','perFraction', ...
    'dvh',[], ...
    'qi',[]);
end

function [doseCubeSample,doseMappingInfo] = mapSampleDoseToReferenceIfNeeded(ct,doseCubeSample,ctScenId,doseMapping)
if doseMapping.enabled
    [doseCubeSample,doseMappingInfo] = matRad_mapDoseToReferenceScenario( ...
        ct,doseCubeSample,ctScenId,doseMapping.refScen);
else
    doseMappingInfo = struct();
    doseMappingInfo.sourceCtScenId = ctScenId;
    doseMappingInfo.referenceCtScenId = doseMapping.refScen;
    doseMappingInfo.method = doseMapping.method;
    doseMappingInfo.mapped = false;
end
end
