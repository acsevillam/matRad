function batches = matRad_makeScenarioDoseRowBatches(rows,batchSize)
% matRad_makeScenarioDoseRowBatches creates row batches with local indices
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

if isempty(rows)
    batches = {};
    return;
end

numBatches = ceil(numel(rows)/batchSize);
batches = cell(numBatches,1);
for b = 1:numBatches
    firstIx = (b - 1)*batchSize + 1;
    lastIx = min(b*batchSize,numel(rows));
    batches{b}.rows = rows(firstIx:lastIx);
    batches{b}.localIx = firstIx:lastIx;
end
end
