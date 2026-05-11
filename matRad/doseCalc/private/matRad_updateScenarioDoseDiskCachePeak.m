function provider = matRad_updateScenarioDoseDiskCachePeak(provider,blockBytes)
% matRad_updateScenarioDoseDiskCachePeak accumulates written cache bytes
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

if ~isfield(provider,'sizeTelemetry') || isempty(provider.sizeTelemetry)
    provider.sizeTelemetry = matRad_initializeScenarioDoseSizeTelemetry();
end

provider.sizeTelemetry.diskCachePeakBytes = ...
    provider.sizeTelemetry.diskCachePeakBytes + double(blockBytes);
end
