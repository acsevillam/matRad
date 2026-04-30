function samplingMemoryEstimate = matRad_estimateSamplingMemory(samplingContext)
% matRad_estimateSamplingMemory estimates sampling output and worker memory
%
% call
%   samplingMemoryEstimate = matRad_estimateSamplingMemory(samplingContext)
%
% input
%   samplingContext:          struct with sampling inputs and diagnostics
%                             needed for memory estimation. Required fields:
%                             ct, stf, cst, cstEval, pln, w, subIx,
%                             samplingCtScenIx, dvhPoints, refGy, refVol,
%                             resultGUInomScen, doseMapping, refScen,
%                             samplingDoseStorageBytes, and numSamples
%
% output
%   samplingMemoryEstimate:   struct with output-storage, raw per-worker,
%                             and component-wise memory estimates in bytes
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

samplingContext = validateSamplingContext(samplingContext);

ct = samplingContext.ct;
stf = samplingContext.stf;
cst = samplingContext.cst;
cstEval = samplingContext.cstEval;
pln = samplingContext.pln;
w = samplingContext.w;
subIx = samplingContext.subIx;
samplingCtScenIx = samplingContext.samplingCtScenIx;
dvhPoints = samplingContext.dvhPoints;
refGy = samplingContext.refGy;
refVol = samplingContext.refVol;
resultGUInomScen = samplingContext.resultGUInomScen;
doseMapping = samplingContext.doseMapping;
refScen = samplingContext.refScen;
samplingDoseStorageBytes = samplingContext.samplingDoseStorageBytes;
numSamples = samplingContext.numSamples;

inputBytes = getVariableBytes(ct) + getVariableBytes(stf) + ...
    getVariableBytes(cst) + getVariableBytes(cstEval) + ...
    getVariableBytes(pln) + getVariableBytes(w) + ...
    getVariableBytes(subIx) + getVariableBytes(samplingCtScenIx) + ...
    getVariableBytes(dvhPoints) + getVariableBytes(refGy) + ...
    getVariableBytes(refVol) + getVariableBytes(doseMapping) + ...
    getVariableBytes(refScen);
sampleDoseBytes = numel(subIx) * 4;
doseResultProxyBytes = getVariableBytes(resultGUInomScen);
sampleResultBytes = estimateSamplingResultBytes(pln,refScen,doseMapping, ...
    resultGUInomScen,numSamples);
doseMappingWorkspaceBytes = 0;
if doseMapping.enabled
    doseMappingWorkspaceBytes = estimateDoseMappingWorkspaceBytes(ct);
end

samplingMemoryEstimate = struct();
samplingMemoryEstimate.estimateBasis = 'nominalScenarioProxy';
samplingMemoryEstimate.numSamples = numSamples;
samplingMemoryEstimate.numVoxels = numel(subIx);
samplingMemoryEstimate.inputBytes = inputBytes;
samplingMemoryEstimate.doseResultProxyBytes = doseResultProxyBytes;
samplingMemoryEstimate.sampleDoseBytes = sampleDoseBytes;
samplingMemoryEstimate.sampleResultBytes = sampleResultBytes;
samplingMemoryEstimate.doseMappingWorkspaceBytes = doseMappingWorkspaceBytes;
samplingMemoryEstimate.rawWorkerBytes = inputBytes + doseResultProxyBytes + ...
    sampleDoseBytes + sampleResultBytes + doseMappingWorkspaceBytes;
samplingMemoryEstimate.sampleDoseStorageBytes = samplingDoseStorageBytes;
samplingMemoryEstimate.sampleResultStorageBytes = sampleResultBytes * numSamples;
samplingMemoryEstimate.mainProcessOutputBytes = samplingMemoryEstimate.sampleDoseStorageBytes + ...
    samplingMemoryEstimate.sampleResultStorageBytes;
samplingMemoryEstimate.workerLimit = [];

end

function samplingContext = validateSamplingContext(samplingContext)
matRad_cfg = MatRad_Config.instance();

if ~isstruct(samplingContext) || ~isscalar(samplingContext)
    matRad_cfg.dispError('samplingContext must be a scalar struct.');
end

requiredFields = {'ct','stf','cst','cstEval','pln','w','subIx', ...
    'samplingCtScenIx','dvhPoints','refGy','refVol','resultGUInomScen', ...
    'doseMapping','refScen','samplingDoseStorageBytes','numSamples'};

for fieldIx = 1:numel(requiredFields)
    fieldName = requiredFields{fieldIx};
    if ~isfield(samplingContext,fieldName)
        matRad_cfg.dispError('samplingContext.%s is required.',fieldName);
    end
end
end

function sampleResultBytes = estimateSamplingResultBytes(pln,refScen,doseMapping,resultGUInomScen,numSamples)
sampleResult = createEmptySamplingResult();
sampleResult.bioParam = pln.bioParam;
sampleResult.ctScen = refScen;
sampleResult.refScen = refScen;
sampleResult.doseMapping = createRepresentativeDoseMappingInfo(refScen,doseMapping);
if numSamples > 0
    representativeScenario = pln.multScen.extractSingleScenario(1);
    sampleResult.relRangeShift = representativeScenario.relRangeShift;
    sampleResult.absRangeShift = representativeScenario.absRangeShift;
    sampleResult.isoShift = representativeScenario.isoShift;
end
if isfield(resultGUInomScen,'dvh')
    sampleResult.dvh = resultGUInomScen.dvh;
end
if isfield(resultGUInomScen,'qi')
    sampleResult.qi = resultGUInomScen.qi;
end
sampleResultBytes = getVariableBytes(sampleResult);
end

function sampleResult = createEmptySamplingResult()
sampleResult = struct( ...
    'bioParam',[], ...
    'ctScen',[], ...
    'refScen',[], ...
    'doseMapping',[], ...
    'relRangeShift',[], ...
    'absRangeShift',[], ...
    'isoShift',[], ...
    'evaluationModeBase','perFraction', ...
    'dvh',[], ...
    'qi',[]);
end

function doseMappingInfo = createRepresentativeDoseMappingInfo(refScen,doseMapping)
doseMappingInfo = struct();
doseMappingInfo.sourceCtScen = refScen;
doseMappingInfo.referenceCtScen = refScen;
doseMappingInfo.method = doseMapping.method;
doseMappingInfo.mapped = false;
if doseMapping.enabled
    doseMappingInfo.mapped = true;
end
end

function bytes = estimateDoseMappingWorkspaceBytes(ct)
% X, Y, Z, DVF components, and mapped dose are allocated as double arrays.
bytes = 7 * prod(double(ct.cubeDim)) * 8;
end

function bytes = getVariableBytes(value) %#ok<INUSD>
info = whos('value');
bytes = info.bytes;
end
