function [targetBatches, oarBatches] = matRad_buildDoseIntervalBatches(ctx, cfg)
% matRad_buildDoseIntervalBatches builds INTERVAL row batches
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

targetBatchSize = ScenarioBatch.Batches.matRad_selectScenarioDoseBatchSize( ...
                                                                           numel(ctx.targetRows), numel(ctx.scenarioDijIx), ctx.numBixels, cfg);
targetBatches = ScenarioBatch.Batches.matRad_buildScenarioDoseRowBatches( ...
                                                                         ctx.targetRows, targetBatchSize);

oarBatchSize = ScenarioBatch.Batches.matRad_selectScenarioDoseBatchSize( ...
                                                                        numel(ctx.oarRows), numel(ctx.scenarioDijIx), ctx.numBixels, cfg);
oarBatches = ScenarioBatch.Batches.matRad_buildScenarioDoseRowBatches(ctx.oarRows, ...
                                                                      oarBatchSize);
end
