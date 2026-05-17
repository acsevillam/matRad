function [expectedRows, expectedBatches, voiRows, voiBatches] = matRad_buildDoseProbBatches(ctx, cfg)
% matRad_buildDoseProbBatches builds probabilistic expected and VOI batches
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

expectedRows = unique([ctx.targetRows(:); ctx.oarRows(:)], 'stable');
expectedBatchSize = ScenarioBatch.Batches.matRad_selectScenarioDoseBatchSize( ...
                                                                             numel(expectedRows), numel(ctx.scenarioDijIx), ctx.numBixels, cfg);
expectedBatches = ScenarioBatch.Batches.matRad_buildScenarioDoseRowBatches( ...
                                                                           expectedRows, expectedBatchSize);

voiRows = matRad_buildDoseProbVoiRows(ctx);
voiBatches = matRad_buildDoseProbVoiRowBatches(voiRows, ctx, cfg);
end

function voiRows = matRad_buildDoseProbVoiRows(ctx)
voiRows = cell(size(ctx.cstDoseGrid, 1), 1);
selectedStructures = [ctx.targetStructures(:); ctx.oarStructures(:)];

for i = 1:numel(selectedStructures)
    rowIx = selectedStructures(i).cstRow;
    voxels = selectedStructures(i).voxelIx(:);
    if isempty(voiRows{rowIx})
        voiRows{rowIx} = voxels;
    else
        voiRows{rowIx} = unique([voiRows{rowIx}; voxels], 'stable');
    end
end
end

function voiBatches = matRad_buildDoseProbVoiRowBatches(voiRows, ctx, cfg)
voiBatches = cell(size(voiRows));
for voiIx = 1:numel(voiRows)
    rows = voiRows{voiIx};
    if isempty(rows)
        voiBatches{voiIx} = {};
        continue
    end

    batchSize = ScenarioBatch.Batches.matRad_selectScenarioDoseBatchSize(numel(rows), ...
                                                                         numel(ctx.scenarioDijIx), ctx.numBixels, cfg);
    voiBatches{voiIx} = ScenarioBatch.Batches.matRad_buildScenarioDoseRowBatches(rows, ...
                                                                                 batchSize);
end
end
