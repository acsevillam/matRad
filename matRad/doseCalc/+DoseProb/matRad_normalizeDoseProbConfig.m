function cfg = matRad_normalizeDoseProbConfig(cfg, matRadCfg)
% matRad_normalizeDoseProbConfig normalizes probabilistic scenario-batch options
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

cfg = ScenarioBatch.Config.matRad_normalizeScenarioDosePrecomputeConfig( ...
                                                                        cfg, matRadCfg, ...
                                                                        ScenarioBatch.Cache.matRad_defaultScenarioDoseCacheRoot('doseProb'));
end
