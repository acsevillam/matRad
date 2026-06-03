function [overheadBytes, details] = matRad_estimateSamplingParallelOverheadBytes( ...
    sampleDoseBytes, sampleResultBytes, forwardDoseWorkerDetails)
% matRad_estimateSamplingParallelOverheadBytes estimates parfor sampling overhead
%
% call
%   overheadBytes = matRad_estimateSamplingParallelOverheadBytes(...)
%   [overheadBytes,details] = matRad_estimateSamplingParallelOverheadBytes(...)
%
% input
%   sampleDoseBytes:           one sampled dose column in bytes
%   sampleResultBytes:         one sampled result proxy in bytes
%   forwardDoseWorkerDetails:  worker estimate component details
%
% output
%   overheadBytes:             extra per-task bytes used for parallel planning
%   details:                   estimate components
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

sampleDoseBytes = matRad_nonnegativeBytes(sampleDoseBytes);
sampleResultBytes = matRad_nonnegativeBytes(sampleResultBytes);

forwardDoseWorkspaceBytes = 0;
if isstruct(forwardDoseWorkerDetails) && ...
        isfield(forwardDoseWorkerDetails, 'forwardDoseWorkspaceBytes')
    forwardDoseWorkspaceBytes = matRad_nonnegativeBytes( ...
        forwardDoseWorkerDetails.forwardDoseWorkspaceBytes);
end

chunkBufferBytesPerTask = sampleDoseBytes;
parforResultCopyBytesPerTask = sampleDoseBytes + sampleResultBytes;
minParallelWorkerOverheadBytes = 2 * 1024^3;
parallelWorkerOverheadRatio = 0.30;
parallelWorkerOverheadBytes = max(minParallelWorkerOverheadBytes, ...
    parallelWorkerOverheadRatio * forwardDoseWorkspaceBytes);

overheadBytes = chunkBufferBytesPerTask + parforResultCopyBytesPerTask + ...
    parallelWorkerOverheadBytes;

details = struct();
details.chunkBufferBytesPerTask = chunkBufferBytesPerTask;
details.parforResultCopyBytesPerTask = parforResultCopyBytesPerTask;
details.parallelWorkerOverheadBytes = parallelWorkerOverheadBytes;
details.parallelWorkerOverheadRatio = parallelWorkerOverheadRatio;
details.minParallelWorkerOverheadBytes = minParallelWorkerOverheadBytes;
details.parallelOverheadBytesPerTask = overheadBytes;
details.forwardDoseWorkspaceBytes = forwardDoseWorkspaceBytes;
end

function value = matRad_nonnegativeBytes(value)
if isempty(value) || ~isnumeric(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 0
    value = 0;
else
    value = double(value);
end
end
