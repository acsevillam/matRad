function streamingSize = matRad_buildScenarioDoseStreamingSize(dij,provider,cfg)
% matRad_buildScenarioDoseStreamingSize summarizes streaming precompute bytes
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

telemetry = matRad_initializeScenarioDoseSizeTelemetry();
if isfield(provider,'sizeTelemetry') && ~isempty(provider.sizeTelemetry)
    telemetry = provider.sizeTelemetry;
end

compactBytes = matRad_variableBytes(dij);
diskCachePeakBytes = double(telemetry.diskCachePeakBytes);
memoryTemporaryPeakBytes = double(telemetry.memoryTemporaryPeakBytes);
if strcmp(cfg.SecondPassStrategy,'disk')
    auxiliaryPeakKind = 'diskCache';
    auxiliaryPeakBytes = diskCachePeakBytes;
else
    auxiliaryPeakKind = 'memoryTemporary';
    auxiliaryPeakBytes = memoryTemporaryPeakBytes;
end

streamingSize = struct();
streamingSize.compactBytes = compactBytes;
streamingSize.auxiliaryPeakBytes = auxiliaryPeakBytes;
streamingSize.totalPrecomputingBytes = compactBytes + auxiliaryPeakBytes;
streamingSize.diskCachePeakBytes = diskCachePeakBytes;
streamingSize.memoryTemporaryPeakBytes = memoryTemporaryPeakBytes;
streamingSize.auxiliaryPeakKind = auxiliaryPeakKind;
streamingSize.secondPassStrategy = cfg.SecondPassStrategy;

dij.streamingSize = streamingSize;
compactBytes = matRad_variableBytes(dij);
streamingSize.compactBytes = compactBytes;
streamingSize.totalPrecomputingBytes = compactBytes + auxiliaryPeakBytes;
end
