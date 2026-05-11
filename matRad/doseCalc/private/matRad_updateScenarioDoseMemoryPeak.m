function provider = matRad_updateScenarioDoseMemoryPeak(provider,varargin)
% matRad_updateScenarioDoseMemoryPeak stores the largest temporary byte estimate
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

bytes = 0;
for i = 1:numel(varargin)
    bytes = bytes + matRad_variableBytes(varargin{i});
end

provider.sizeTelemetry.memoryTemporaryPeakBytes = max( ...
    provider.sizeTelemetry.memoryTemporaryPeakBytes,double(bytes));
end
