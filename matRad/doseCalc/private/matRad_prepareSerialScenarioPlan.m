function pln_s = matRad_prepareSerialScenarioPlan(pln,scenarioId)
% matRad_prepareSerialScenarioPlan prepares one serial scenario plan
%
% call
%   pln_s = matRad_prepareSerialScenarioPlan(pln,scenarioId)
%
% input
%   pln:        matRad plan struct with pln.multScen as matRad_ScenarioModel
%   scenarioId: public scenario id selected from pln.multScen.scenarioIds()
%
% output
%   pln_s:      single-scenario worker-safe plan with nested parallelism off
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

pln_s = matRad_selectSingleScenarioPlan(pln,scenarioId);
pln_s = matRad_makeWorkerSafePlan(pln_s);
end
