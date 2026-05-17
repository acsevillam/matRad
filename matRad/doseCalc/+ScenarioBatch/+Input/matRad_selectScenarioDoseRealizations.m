function scenarioInfo = matRad_selectScenarioDoseRealizations(ct, pln, cfg, matRadCfg)
% matRad_selectScenarioDoseRealizations selects active scenario-batch scenarios
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

if ~isfield(pln, 'multScen') || isempty(pln.multScen) || ...
        ~isa(pln.multScen, 'matRad_ScenarioModel')
    matRadCfg.dispError('Scenario dose precomputation requires pln.multScen.');
end

refScen = 1;
if isfield(ct, 'refScen') && ~isempty(ct.refScen)
    refScen = ct.refScen;
end
if isfield(cfg, 'refScen') && ~isempty(cfg.refScen)
    refScen = cfg.refScen;
end
matRad_validatePositiveInteger(refScen, 'refScen', matRadCfg);

scenarioIds = pln.multScen.scenarioIds();
scenarioIds = scenarioIds(:);
scenarioDijIx = arrayfun(@(id) pln.multScen.getDijScenarioIndex(id), scenarioIds);
scenarioDijIx = scenarioDijIx(:);
scenarioCtScenIds = arrayfun(@(id) pln.multScen.getCtScenario(id), scenarioIds);
scenarioCtScenIds = scenarioCtScenIds(:);
selectedCtScenIds = matRad_selectScen4DCtScenarios(pln, scenarioCtScenIds, matRadCfg);
scenarioMask = ismember(scenarioCtScenIds, selectedCtScenIds);
scenarioMask = scenarioMask(:);

if ~any(scenarioMask)
    matRadCfg.dispError('pln.propOpt.scen4D does not select any active CT scenarios.');
end

scenarioInfo.scenarioIds = scenarioIds(scenarioMask);
scenarioInfo.scenarioDijIx = scenarioDijIx(scenarioMask);
scenarioInfo.scenarioCtScenIds = scenarioCtScenIds(scenarioMask);
scenarioInfo.scenarioWeights = matRad_buildSelectedScenarioWeights(pln.multScen, scenarioMask, matRadCfg);
scenarioInfo.scenarioMask = scenarioMask;
scenarioInfo.refScen = refScen;
end

function selectedCtScenIds = matRad_selectScen4DCtScenarios(pln, scenarioCtScenIds, matRadCfg)
activeCtScenIds = unique(scenarioCtScenIds(:), 'stable');
selectedCtScenIds = 1;

if ~isfield(pln, 'propOpt') || ~isstruct(pln.propOpt) || ...
        ~isfield(pln.propOpt, 'scen4D') || isempty(pln.propOpt.scen4D)
    matRad_validateSelectedCtScenIds(selectedCtScenIds, activeCtScenIds, matRadCfg);
    return
end

scen4D = pln.propOpt.scen4D;
if isstring(scen4D) && isscalar(scen4D)
    scen4D = char(scen4D);
end

if ischar(scen4D)
    scen4D = strtrim(scen4D);
    if strcmpi(scen4D, 'all')
        selectedCtScenIds = activeCtScenIds;
        return
    end
    matRadCfg.dispError(['pln.propOpt.scen4D must be ''all'' or a numeric ', ...
                         'vector of positive integer CT scenario ids.']);
end

if ~isnumeric(scen4D) || isempty(scen4D) || any(~isfinite(scen4D(:))) || ...
        any(scen4D(:) < 1) || any(fix(scen4D(:)) ~= scen4D(:))
    matRadCfg.dispError(['pln.propOpt.scen4D must be ''all'' or a numeric ', ...
                         'vector of positive integer CT scenario ids.']);
end

selectedCtScenIds = unique(double(scen4D(:)), 'stable');
matRad_validateSelectedCtScenIds(selectedCtScenIds, activeCtScenIds, matRadCfg);
end

function matRad_validateSelectedCtScenIds(selectedCtScenIds, activeCtScenIds, matRadCfg)
missingCtScenIds = selectedCtScenIds(~ismember(selectedCtScenIds, activeCtScenIds));
if ~isempty(missingCtScenIds)
    matRadCfg.dispError('Scenario dose scen4D selection includes inactive CT scenario ids: %s.', ...
                        mat2str(missingCtScenIds(:)'));
end
end

function scenarioWeights = matRad_buildSelectedScenarioWeights(multScen, scenarioMask, matRadCfg)
if matRad_hasFieldOrProp(multScen, 'scenWeight') && ~isempty(multScen.scenWeight)
    rawWeights = multScen.scenWeight(:);
elseif matRad_hasFieldOrProp(multScen, 'scenProb') && ~isempty(multScen.scenProb)
    rawWeights = multScen.scenProb(:);
else
    matRadCfg.dispError('Scenario model must provide scenWeight or scenProb.');
end

if numel(rawWeights) ~= numel(scenarioMask)
    matRadCfg.dispError('Number of scenario weights is inconsistent with active scenarios.');
end
scenarioWeights = rawWeights(scenarioMask);

scenarioWeights = scenarioWeights(:);
scenarioWeights(~isfinite(scenarioWeights) | scenarioWeights < 0) = 0;
if sum(scenarioWeights) <= 0
    matRadCfg.dispError('Scenario weights must contain at least one positive finite value.');
end
scenarioWeights = scenarioWeights ./ sum(scenarioWeights);
end

function tf = matRad_hasFieldOrProp(value, fieldName)
tf = (isobject(value) && isprop(value, fieldName)) || ...
    (isstruct(value) && isfield(value, fieldName));
end

function value = matRad_validatePositiveInteger(value, name, matRadCfg)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value < 1 || fix(value) ~= value
    matRadCfg.dispError('%s must be a positive integer scalar.', name);
end
value = double(value);
end
