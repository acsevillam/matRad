function pln_s = matRad_selectSingleScenarioPlan(pln,scenarioId)
% matRad_selectSingleScenarioPlan returns a plan with one active scenario
%
% call
%   pln_s = matRad_selectSingleScenarioPlan(pln,scenarioId)
%
% input
%   pln:        matRad plan struct with pln.multScen as matRad_ScenarioModel
%   scenarioId: public scenario id selected from pln.multScen.scenarioIds()
%
% output
%   pln_s:      copy of pln whose multScen contains only scenarioId
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

matRad_cfg = MatRad_Config.instance();

if ~isfield(pln,'multScen') || isempty(pln.multScen) || ...
        ~isa(pln.multScen,'matRad_ScenarioModel')
    matRad_cfg.dispError('Single-scenario dose calculation requires pln.multScen as matRad_ScenarioModel.');
end

pln_s = pln;
pln_s.multScen = pln.multScen.extractSingleScenario(scenarioId);

if isfield(pln_s,'propOpt') && isstruct(pln_s.propOpt) && ...
        isfield(pln_s.propOpt,'scen4D')
    pln_s.propOpt = rmfield(pln_s.propOpt,'scen4D');
end

end
