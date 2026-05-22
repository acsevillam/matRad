function [workerBytes, details] = matRad_estimateSamplingForwardDoseWorkerBytes( ...
    samplingContext, inputBytes, doseCubeBytes, doseMappingWorkspaceBytes, ...
    minForwardDoseWorkerMemoryBytes)
% matRad_estimateSamplingForwardDoseWorkerBytes estimates sampling worker memory
%
% call
%   workerBytes = matRad_estimateSamplingForwardDoseWorkerBytes(...)
%   [workerBytes,details] = matRad_estimateSamplingForwardDoseWorkerBytes(...)
%
% input
%   samplingContext: sampling context used by one worker
%   inputBytes:      serialized input proxy size
%   doseCubeBytes:   one dose cube proxy size
%
% output
%   workerBytes:     conservative worker memory estimate in bytes
%   details:         estimate components
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

if nargin < 5 || isempty(minForwardDoseWorkerMemoryBytes)
    minForwardDoseWorkerMemoryBytes = 16 * 1024^3;
end

numBeams = matRad_samplingBeamCount(samplingContext);
ctVoxelCount = prod(double(samplingContext.ct.cubeDim));
selectedVoxelCount = numel(samplingContext.subIx);

inputBytes = matRad_nonnegativeBytes(inputBytes);
doseCubeBytes = matRad_nonnegativeBytes(doseCubeBytes);
doseMappingWorkspaceBytes = matRad_nonnegativeBytes(doseMappingWorkspaceBytes);
minForwardDoseWorkerMemoryBytes = matRad_nonnegativeBytes(minForwardDoseWorkerMemoryBytes);

beamWorkspaceBytes = doseCubeBytes * max(1, numBeams) * 2;
voxelWorkspaceBytes = max(ctVoxelCount, selectedVoxelCount) * 8 * 4;
matlabWorkerRuntimeBytes = 2 * 1024^3;

forwardDoseWorkspaceBytes = max([minForwardDoseWorkerMemoryBytes, ...
                                 beamWorkspaceBytes, ...
                                 voxelWorkspaceBytes, ...
                                 matlabWorkerRuntimeBytes]);

workerBytes = inputBytes + doseCubeBytes + doseMappingWorkspaceBytes + ...
    forwardDoseWorkspaceBytes;

details = struct();
details.estimateBasis = 'samplingForwardDose';
details.inputBytes = inputBytes;
details.doseCubeBytes = doseCubeBytes;
details.doseMappingWorkspaceBytes = doseMappingWorkspaceBytes;
details.forwardDoseWorkspaceBytes = forwardDoseWorkspaceBytes;
details.minForwardDoseWorkerMemoryBytes = minForwardDoseWorkerMemoryBytes;
details.beamWorkspaceBytes = beamWorkspaceBytes;
details.voxelWorkspaceBytes = voxelWorkspaceBytes;
details.matlabWorkerRuntimeBytes = matlabWorkerRuntimeBytes;
details.numBeams = numBeams;
details.ctVoxelCount = ctVoxelCount;
details.selectedVoxelCount = selectedVoxelCount;
end

function numBeams = matRad_samplingBeamCount(samplingContext)
numBeams = 1;
if isstruct(samplingContext) && isfield(samplingContext, 'stf') && ...
        ~isempty(samplingContext.stf)
    numBeams = numel(samplingContext.stf);
end
end

function value = matRad_nonnegativeBytes(value)
if isempty(value) || ~isnumeric(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 0
    value = 0;
else
    value = double(value);
end
end
