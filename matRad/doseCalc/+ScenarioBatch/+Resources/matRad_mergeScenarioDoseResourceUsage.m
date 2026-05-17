function provider = matRad_mergeScenarioDoseResourceUsage(provider, resourceBlocks)
% matRad_mergeScenarioDoseResourceUsage merges worker scenario-batch resource usage
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

if ~isfield(provider, 'resourceUsage') || isempty(provider.resourceUsage)
    provider.resourceUsage = ScenarioBatch.Resources.matRad_buildScenarioDoseResourceUsage();
end

for i = 1:numel(resourceBlocks)
    resourceUsage = resourceBlocks{i};
    if isstruct(resourceUsage) && isfield(resourceUsage, 'resourceUsage')
        resourceUsage = resourceUsage.resourceUsage;
    end
    if isempty(resourceUsage)
        continue
    end

    if isfield(resourceUsage, 'diskCachePeakBytes')
        provider.resourceUsage.diskCachePeakBytes = ...
            provider.resourceUsage.diskCachePeakBytes + ...
            double(resourceUsage.diskCachePeakBytes);
    end

    if isfield(resourceUsage, 'memoryTemporaryPeakBytes')
        provider.resourceUsage.memoryTemporaryPeakBytes = max( ...
                                                              provider.resourceUsage.memoryTemporaryPeakBytes, ...
                                                              double(resourceUsage.memoryTemporaryPeakBytes));
    end
end
end
