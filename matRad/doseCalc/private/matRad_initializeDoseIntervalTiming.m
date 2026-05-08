function timing = matRad_initializeDoseIntervalTiming(kind,varargin)
% matRad_initializeDoseIntervalTiming creates dose interval timing diagnostics
%
% call
%   timing = matRad_initializeDoseIntervalTiming('interval',intervalMode,cfg, ...
%       numTargetVoxels,numOarVoxels,numScenarios,numBixels)
%   timing = matRad_initializeDoseIntervalTiming('stage',stageName,numVoxels, ...
%       batchSize,numBatches)
%   timing = matRad_initializeDoseIntervalTiming('batch',numVoxels)
%
% input
%   kind:          timing struct kind, one of 'interval', 'stage', or 'batch'
%   varargin:      kind-specific initialization values
%
% output
%   timing:        scalar struct with zero-initialized timing diagnostics
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

kind = normalizeTimingKind(kind);
switch kind
    case 'interval'
        timing = initializeIntervalTiming(varargin{:});
    case 'stage'
        timing = initializeStageTiming(varargin{:});
    case 'batch'
        timing = initializeBatchTiming(varargin{:});
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Unsupported dose interval timing kind "%s".',kind);
end
end

function kind = normalizeTimingKind(kind)
if isstring(kind) && isscalar(kind)
    kind = char(kind);
end
if ~ischar(kind) || isempty(kind)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Dose interval timing kind must be non-empty text.');
end
kind = lower(strtrim(kind));
end

function timing = initializeIntervalTiming(intervalMode,cfg,numTargetVoxels, ...
    numOarVoxels,numScenarios,numBixels)
timing = struct();
timing.intervalMode = intervalMode;
timing.useParallelRequested = logical(cfg.UseParallel);
timing.numTargetVoxels = numTargetVoxels;
timing.numOarVoxels = numOarVoxels;
timing.numScenarios = numScenarios;
timing.numBixels = numBixels;
timing.totalSeconds = 0;
end

function timing = initializeStageTiming(stageName,numVoxels,batchSize,numBatches)
timing = struct();
timing.stage = stageName;
timing.numVoxels = numVoxels;
timing.batchSize = batchSize;
timing.numBatches = numBatches;
timing.extractMapSeconds = 0;
timing.centerAccumSeconds = 0;
timing.radiusMultiplySeconds = 0;
timing.factorSeconds = 0;
timing.centeredRowsSeconds = 0;
timing.parallelSetupSeconds = 0;
timing.parallelComputeWallSeconds = 0;
timing.serialAssemblySeconds = 0;
timing.wallSeconds = 0;
timing.parallelEnabled = false;
end

function timing = initializeBatchTiming(numVoxels)
timing = initializeStageTiming('batch',numVoxels,numVoxels,1);
end
