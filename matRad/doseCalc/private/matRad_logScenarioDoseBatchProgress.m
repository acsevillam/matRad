function matRad_logScenarioDoseBatchProgress(matRad_cfg,cfg,stageName,batchIx,numBatches)
% matRad_logScenarioDoseBatchProgress logs streaming batch progress
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

if numBatches == 0
    return;
end

progressStep = max(1,ceil(numBatches/10));
if isfield(cfg,'ProgressLevel') && ...
        (strcmp(cfg.ProgressLevel,'detailed') || batchIx == 1 || ...
        batchIx == numBatches || mod(batchIx,progressStep) == 0)
    matRad_cfg.dispInfo('matRad: %s progress %d/%d (%.0f%%).\n', ...
        stageName,batchIx,numBatches,100*batchIx/numBatches);
end
end
