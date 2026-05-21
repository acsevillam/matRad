function workerMemoryBytes = matRad_estimateScenarioDoseWorkerMemoryBytes( ...
                                                                          parallelProvider, originalProvider, ctx)
% matRad_estimateScenarioDoseWorkerMemoryBytes estimates worker memory
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

providerBytes = ScenarioBatch.Resources.matRad_estimateVariableBytes(parallelProvider);
preloadedDijBytes = 0;
if isfield(originalProvider, 'preloadedDij') && ~isempty(originalProvider.preloadedDij)
    preloadedDijBytes = ScenarioBatch.Resources.matRad_estimateVariableBytes( ...
                                                                             originalProvider.preloadedDij);
end

numVoxels = max(1, double(ctx.numVoxels));
numBixels = max(1, double(ctx.numBixels));
matrixFloorBytes = 0.25 * numVoxels * numBixels * 8;
minScenarioWorksetBytes = 4 * 1024^3;
scenarioDijBytes = max([preloadedDijBytes, matrixFloorBytes, ...
                        minScenarioWorksetBytes]);

numRows = max([1; numel(ctx.targetRows); numel(ctx.oarRows)]);
sparseWorkspaceFillFactor = 0.05;
rowWorkspaceBytes = ScenarioBatch.Resources.matRad_estimateSparseMatrixBytes( ...
                                                                             numRows, numBixels, sparseWorkspaceFillFactor);
temporaryWorkspaceBytes = max(128 * 1024^2, double(numRows) * 8 * 10);
workerMemoryBytes = providerBytes + scenarioDijBytes + rowWorkspaceBytes + ...
    temporaryWorkspaceBytes;
end
