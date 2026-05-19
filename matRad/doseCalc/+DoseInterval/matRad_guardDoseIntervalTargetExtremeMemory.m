function useParallel = matRad_guardDoseIntervalTargetExtremeMemory(ctx, targetBatches, cfg, useParallel, matRadCfg)
% matRad_guardDoseIntervalTargetExtremeMemory checks target extreme memory
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
estimatedMB = matRad_estimateIntervalTargetExtremeMemoryMB(ctx, targetBatches, ...
                                                           useParallel);

matRadCfg.dispInfo(['matRad: Estimated INTERVAL target extreme radius memory ', ...
                    'is %.2f MB (limit %.2f MB).\n'], ...
                   estimatedMB, memoryLimitMB);

if estimatedMB <= memoryLimitMB
    return
end

if useParallel
    serialEstimatedMB = matRad_estimateIntervalTargetExtremeMemoryMB(ctx, ...
                                                                     targetBatches, false);
    matRadCfg.dispInfo(['matRad: Estimated serial INTERVAL target extreme ', ...
                        'radius memory is %.2f MB (limit %.2f MB).\n'], ...
                       serialEstimatedMB, memoryLimitMB);

    if serialEstimatedMB <= memoryLimitMB
        matRadCfg.dispWarning(['Parallel INTERVAL target extreme radius estimated ', ...
                               'memory is %.2f MB, which exceeds MemoryLimitMB %.2f MB. ', ...
                               'Falling back to serial target extreme radius for ', ...
                               'the second pass.\n'], ...
                              estimatedMB, memoryLimitMB);
        useParallel = false;
        return
    end

    matRadCfg.dispError(['INTERVAL target extreme radius estimated memory is %.2f MB ', ...
                         'in parallel and %.2f MB in serial, which exceeds MemoryLimitMB ', ...
                         '%.2f MB. Increase cfg.MemoryLimitMB, reduce the number of ', ...
                         'bixels, scenarios, target voxels, or batch rows.'], ...
                        estimatedMB, serialEstimatedMB, memoryLimitMB);
end

matRadCfg.dispError(['INTERVAL target extreme radius estimated memory is %.2f MB, ', ...
                     'which exceeds MemoryLimitMB %.2f MB. Increase cfg.MemoryLimitMB, ', ...
                     'reduce the number of bixels, scenarios, target voxels, or batch rows.'], ...
                    estimatedMB, memoryLimitMB);
end

function estimatedMB = matRad_estimateIntervalTargetExtremeMemoryMB(ctx, targetBatches, useParallel)
bytesPerDouble = 8;
numScenarios = numel(ctx.scenarioDijIx);
numTargetRows = numel(ctx.targetRows);
numBixels = ctx.numBixels;
maxBatchRows = matRad_getMaxBatchRows(targetBatches);

if numScenarios == 0 || numTargetRows == 0 || numBixels == 0 || maxBatchRows == 0
    estimatedMB = 0;
    return
end

% Target extreme stores sparse row-wise deltas. Use the same conservative
% fill heuristic used for scenario-dose worker row workspaces.
sparseWorkspaceFillFactor = 0.05;
targetDeltaBytes = numTargetRows * numBixels * bytesPerDouble * ...
    sparseWorkspaceFillFactor;
batchWorkBytes = 2 * maxBatchRows * numBixels * bytesPerDouble * ...
    sparseWorkspaceFillFactor;

parallelBytes = 0;
if useParallel
    centerRowsBytes = targetDeltaBytes;
    parallelResultBytes = numScenarios * targetDeltaBytes;
    parallelBytes = centerRowsBytes + parallelResultBytes;
end

estimatedMB = (targetDeltaBytes + batchWorkBytes + parallelBytes) / 1e6;
end

function maxBatchRows = matRad_getMaxBatchRows(batches)
maxBatchRows = 0;
for batchIx = 1:numel(batches)
    maxBatchRows = max(maxBatchRows, numel(batches{batchIx}.rows));
end
end
