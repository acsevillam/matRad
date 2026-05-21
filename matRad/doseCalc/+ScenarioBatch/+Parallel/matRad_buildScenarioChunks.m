function chunks = matRad_buildScenarioChunks(numScenarios, chunkSize)
% matRad_buildScenarioChunks partitions scenario indices into chunks
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

matRad_validatePositiveInteger(numScenarios, 'numScenarios');
matRad_validatePositiveInteger(chunkSize, 'chunkSize');

chunkSize = min(double(chunkSize), double(numScenarios));
numChunks = ceil(double(numScenarios) / chunkSize);
chunks = cell(numChunks, 1);

for chunkIx = 1:numChunks
    firstScenario = (chunkIx - 1) * chunkSize + 1;
    lastScenario = min(double(numScenarios), firstScenario + chunkSize - 1);
    chunks{chunkIx} = firstScenario:lastScenario;
end
end

function matRad_validatePositiveInteger(value, valueName)
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
     round(value) == value && value >= 1)
    matRadCfg = MatRad_Config.instance();
    matRadCfg.dispError('%s must be a positive integer scalar.', valueName);
end
end
