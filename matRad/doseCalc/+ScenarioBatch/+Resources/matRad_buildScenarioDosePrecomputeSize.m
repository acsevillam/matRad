function precomputeSize = matRad_buildScenarioDosePrecomputeSize(dij, provider, cfg)
% matRad_buildScenarioDosePrecomputeSize summarizes scenario-batch precompute bytes
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

resourceUsage = ScenarioBatch.Resources.matRad_buildScenarioDoseResourceUsage();
if isfield(provider, 'resourceUsage') && ~isempty(provider.resourceUsage)
    resourceUsage = provider.resourceUsage;
end

compactBytes = ScenarioBatch.Resources.matRad_estimateVariableBytes(dij);
diskCachePeakBytes = double(resourceUsage.diskCachePeakBytes);
memoryTemporaryPeakBytes = double(resourceUsage.memoryTemporaryPeakBytes);
if strcmp(cfg.SecondPassStrategy, 'disk')
    auxiliaryPeakKind = 'diskCache';
    auxiliaryPeakBytes = diskCachePeakBytes;
else
    auxiliaryPeakKind = 'memoryTemporary';
    auxiliaryPeakBytes = memoryTemporaryPeakBytes;
end

precomputeSize = struct();
precomputeSize.compactBytes = compactBytes;
precomputeSize.auxiliaryPeakBytes = auxiliaryPeakBytes;
precomputeSize.totalPrecomputingBytes = compactBytes + auxiliaryPeakBytes;
precomputeSize.diskCachePeakBytes = diskCachePeakBytes;
precomputeSize.memoryTemporaryPeakBytes = memoryTemporaryPeakBytes;
precomputeSize.auxiliaryPeakKind = auxiliaryPeakKind;
precomputeSize.secondPassStrategy = cfg.SecondPassStrategy;

dij.precomputeSize = precomputeSize;
compactBytes = ScenarioBatch.Resources.matRad_estimateVariableBytes(dij);
precomputeSize.compactBytes = compactBytes;
precomputeSize.totalPrecomputingBytes = compactBytes + auxiliaryPeakBytes;
end
