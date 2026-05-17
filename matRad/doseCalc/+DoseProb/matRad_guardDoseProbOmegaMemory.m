function matRad_guardDoseProbOmegaMemory(ctx, voiBatches, cfg, useParallel, matRadCfg)
% matRad_guardDoseProbOmegaMemory checks probabilistic omega memory
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

memoryLimitMB = cfg.MemoryLimitMB;
estimatedMB = matRad_estimateDoseProbOmegaMemoryMB(ctx, voiBatches, cfg, ...
                                                   useParallel);

matRadCfg.dispInfo(['matRad: Estimated probabilistic omega memory is %.2f MB ', ...
                    '(limit %.2f MB).\n'], estimatedMB, memoryLimitMB);

if estimatedMB <= memoryLimitMB
    return
end

if useParallel
    matRadCfg.dispError(['Probabilistic omega estimated memory is %.2f MB, which ', ...
                         'exceeds MemoryLimitMB %.2f MB. Increase cfg.MemoryLimitMB, ', ...
                         'reduce the number of bixels, scenarios, VOIs, or batch rows, ', ...
                         'or disable parallel probabilistic omega calculation.'], ...
                        estimatedMB, memoryLimitMB);
end

matRadCfg.dispError(['Probabilistic omega estimated memory is %.2f MB, which ', ...
                     'exceeds MemoryLimitMB %.2f MB. Increase cfg.MemoryLimitMB, ', ...
                     'reduce the number of bixels, scenarios, VOIs, or batch rows.'], ...
                    estimatedMB, memoryLimitMB);
end

function estimatedMB = matRad_estimateDoseProbOmegaMemoryMB(ctx, voiBatches, ...
                                                            cfg, useParallel)
bytesPerDouble = 8;
numScenarios = numel(ctx.scenarioDijIx);
numBixels = ctx.numBixels;
[numActiveVois, maxBatchRows] = matRad_getDoseProbOmegaBatchStats(voiBatches);

if numActiveVois == 0 || maxBatchRows == 0
    estimatedMB = 0;
    return
end

omegaMatrixBytes = numBixels^2 * bytesPerDouble;
omegaOutputBytes = numActiveVois * omegaMatrixBytes;
omegaWorkBytes = 2 * omegaMatrixBytes;

if strcmp(cfg.SecondPassStrategy, 'disk') && ~useParallel
    centeredRowsBytes = numScenarios * maxBatchRows * numBixels * bytesPerDouble;
else
    centeredRowsBytes = maxBatchRows * numBixels * bytesPerDouble;
end

parallelResultBytes = 0;
if useParallel
    parallelResultBytes = numScenarios * numActiveVois * omegaMatrixBytes;
end

estimatedMB = (omegaOutputBytes + omegaWorkBytes + centeredRowsBytes + ...
               parallelResultBytes) / 1e6;
end

function [numActiveVois, maxBatchRows] = matRad_getDoseProbOmegaBatchStats(voiBatches)
numActiveVois = 0;
maxBatchRows = 0;

for voiIx = 1:numel(voiBatches)
    batches = voiBatches{voiIx};
    if isempty(batches)
        continue
    end

    numActiveVois = numActiveVois + 1;
    for batchIx = 1:numel(batches)
        maxBatchRows = max(maxBatchRows, numel(batches{batchIx}.rows));
    end
end
end
