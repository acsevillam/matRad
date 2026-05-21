function bytes = matRad_estimateSparseMatrixBytes(numRows, numCols, fillFactor)
% matRad_estimateSparseMatrixBytes estimates sparse double matrix storage
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

if nargin < 3 || isempty(fillFactor)
    fillFactor = 0.05;
end

numRows = max(0, double(numRows));
numCols = max(0, double(numCols));
fillFactor = min(1, max(0, double(fillFactor)));

nnzEstimate = max(1, ceil(numRows * numCols * fillFactor));
bytesPerDouble = 8;
bytesPerIndex = 8;

bytes = nnzEstimate * (bytesPerDouble + bytesPerIndex) + ...
    (numCols + 1) * bytesPerIndex;
end
