function plnS = matRad_buildSingleScenarioPlan(pln, scenarioId)
% matRad_buildSingleScenarioPlan builds one scenario realization plan
%
% call
%   plnS = ScenarioBatch.Scenarios.matRad_buildSingleScenarioPlan(pln,scenarioId)
%
% input
%   pln:        matRad plan struct with pln.multScen as matRad_ScenarioModel
%   scenarioId: public scenario id selected from pln.multScen.scenarioIds()
%
% output
%   plnS:      worker-safe plan for one active scenario realization
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

matRadCfg = MatRad_Config.instance();
if ~isfield(pln, 'multScen') || isempty(pln.multScen) || ...
        ~isa(pln.multScen, 'matRad_ScenarioModel')
    matRadCfg.dispError('Single-scenario dose calculation requires pln.multScen as matRad_ScenarioModel.');
end

plnS = pln;
plnS.multScen = pln.multScen.extractSingleScenario(scenarioId);
if isfield(plnS, 'propOpt') && isstruct(plnS.propOpt) && ...
        isfield(plnS.propOpt, 'scen4D')
    plnS.propOpt = rmfield(plnS.propOpt, 'scen4D');
end

plnS = ScenarioBatch.Worker.matRad_sanitizeWorkerPlan(plnS);
end
