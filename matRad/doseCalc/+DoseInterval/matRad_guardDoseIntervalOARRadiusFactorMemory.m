function matRad_guardDoseIntervalOARRadiusFactorMemory(numScenarios, numBixels, cfg, matRadCfg)
% matRad_guardDoseIntervalOARRadiusFactorMemory checks INTERVAL3 memory
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
estimatedMB = matRad_estimateIntervalOARRadiusFactorMemoryMB(numScenarios, ...
                                                             numBixels, cfg);

matRadCfg.dispInfo(['matRad: Estimated INTERVAL3 OAR radius factor memory ', ...
                    'per voxel is %.2f MB (limit %.2f MB).\n'], estimatedMB, memoryLimitMB);

if estimatedMB > memoryLimitMB
    matRadCfg.dispError(['INTERVAL3 OAR radius factor estimated memory per ', ...
                         'voxel is %.2f MB, which exceeds MemoryLimitMB %.2f MB. Increase ', ...
                         'cfg.MemoryLimitMB, reduce cfg.KMax, reduce the number of scenarios ', ...
                         'or bixels, or use INTERVAL2.'], estimatedMB, memoryLimitMB);
end
end

function estimatedMB = matRad_estimateIntervalOARRadiusFactorMemoryMB(numScenarios, numBixels, cfg)
bytesPerDouble = 8;
kMax = min([cfg.KMax, numBixels, numScenarios]);

scenarioMatrixBytes = 3 * numScenarios * numBixels * bytesPerDouble;
scenarioGramBytes = 3 * numScenarios^2 * bytesPerDouble;
factorBytes = (3 * numBixels * kMax + kMax^2) * bytesPerDouble;

estimatedMB = (scenarioMatrixBytes + scenarioGramBytes + ...
               factorBytes) / 1e6;
end
