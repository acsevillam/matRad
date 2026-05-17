function localProvider = matRad_buildScenarioDoseWorkerProvider(provider)
% matRad_buildScenarioDoseWorkerProvider builds a worker-local provider copy
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

localProvider = provider;
if ~isfield(localProvider, 'resourceUsage') || isempty(localProvider.resourceUsage)
    localProvider.resourceUsage = ScenarioBatch.Resources.matRad_buildScenarioDoseResourceUsage();
end
end
