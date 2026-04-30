function matRad_logEstimatedSamplingMemory(samplingMemoryEstimate,memoryEstimate,matRad_cfg)
% matRad_logEstimatedSamplingMemory logs sampling memory diagnostics
%
% call
%   matRad_logEstimatedSamplingMemory(samplingMemoryEstimate,memoryEstimate,matRad_cfg)
%
% input
%   samplingMemoryEstimate: sampling output and worker memory estimate
%   memoryEstimate:         memory-limited worker estimate, or empty
%   matRad_cfg:             MatRad_Config instance used for logging
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

matRad_cfg.dispInfo(['matRad: EstimatedSamplingMemory: ', ...
    formatMemorySize(samplingMemoryEstimate.mainProcessOutputBytes), ...
    ' output storage (dose matrix ', ...
    formatMemorySize(samplingMemoryEstimate.sampleDoseStorageBytes), ...
    ', result structs ', ...
    formatMemorySize(samplingMemoryEstimate.sampleResultStorageBytes), ...
    '), ', formatMemorySize(samplingMemoryEstimate.rawWorkerBytes), ...
    ' raw per active worker/sample.\n']);

matRad_cfg.dispInfo(['matRad: EstimatedSamplingMemory worker components: inputs ', ...
    formatMemorySize(samplingMemoryEstimate.inputBytes), ...
    ', dose result proxy ', ...
    formatMemorySize(samplingMemoryEstimate.doseResultProxyBytes), ...
    ', sampled dose column ', ...
    formatMemorySize(samplingMemoryEstimate.sampleDoseBytes), ...
    ', sample result ', ...
    formatMemorySize(samplingMemoryEstimate.sampleResultBytes), ...
    ', dose mapping workspace ', ...
    formatMemorySize(samplingMemoryEstimate.doseMappingWorkspaceBytes), '.\n']);

if isempty(memoryEstimate) || ~isstruct(memoryEstimate) || isempty(memoryEstimate.maxWorkers)
    return;
end

memorySource = '';
if isfield(memoryEstimate,'source') && ~isempty(memoryEstimate.source)
    memorySource = [' (', memoryEstimate.source, ')'];
end

workerFloorBytes = [];
if isfield(memoryEstimate,'minWorkerMemoryBytes')
    workerFloorBytes = memoryEstimate.minWorkerMemoryBytes;
end
safetyFactor = 1;
if isfield(memoryEstimate,'safetyFactor')
    safetyFactor = memoryEstimate.safetyFactor;
end

concurrentWorkerBytes = memoryEstimate.workerBytes * memoryEstimate.maxWorkers;
matRad_cfg.dispInfo(['matRad: EstimatedSamplingMemory worker limit: ', ...
    formatMemorySize(memoryEstimate.usableBytes), ' usable', memorySource, ...
    ', max(raw ', formatMemorySize(memoryEstimate.rawWorkerBytes), ...
    ', floor ', formatMemorySize(workerFloorBytes), ') * safety ', ...
    num2str(safetyFactor,'%.2f'), ' = ', ...
    formatMemorySize(memoryEstimate.workerBytes), ...
    ' effective per worker, max ', num2str(memoryEstimate.maxWorkers), ...
    ' worker(s), ', formatMemorySize(concurrentWorkerBytes), ...
    ' estimated concurrent worker memory.\n']);

end

function text = formatMemorySize(bytes)
if isempty(bytes) || ~isnumeric(bytes) || ~isscalar(bytes) || ~isfinite(bytes)
    text = 'n/a';
else
    text = [num2str(bytes / 1024^3,'%.2f'), ' GB'];
end
end
