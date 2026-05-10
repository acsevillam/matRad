function scenarioInfo = matRad_resolveStreamingScenarioSelection(ct,pln,cfg,matRad_cfg)
% matRad_resolveStreamingScenarioSelection selects active streaming scenarios
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

if ~isfield(pln,'multScen') || isempty(pln.multScen) || ...
        ~isa(pln.multScen,'matRad_ScenarioModel')
    matRad_cfg.dispError('Streaming scenario dose calculation requires pln.multScen.');
end

refScen = 1;
if isfield(ct,'refScen') && ~isempty(ct.refScen)
    refScen = ct.refScen;
end
if isfield(cfg,'refScen') && ~isempty(cfg.refScen)
    refScen = cfg.refScen;
end
validatePositiveInteger(refScen,'refScen',matRad_cfg);

scenarioIds = pln.multScen.scenarioIds();
scenarioDijIx = arrayfun(@(id) pln.multScen.getDijScenarioIndex(id),scenarioIds);
scenarioCtScenIds = arrayfun(@(id) pln.multScen.getCtScenario(id),scenarioIds);
selectedCtScenIds = resolveScen4DSelection(pln,scenarioCtScenIds,matRad_cfg);
scenarioMask = ismember(scenarioCtScenIds,selectedCtScenIds);

if ~any(scenarioMask)
    matRad_cfg.dispError('pln.propOpt.scen4D does not select any active CT scenarios.');
end

scenarioInfo.scenarioIds = scenarioIds(scenarioMask);
scenarioInfo.scenarioDijIx = scenarioDijIx(scenarioMask);
scenarioInfo.scenarioCtScenIds = scenarioCtScenIds(scenarioMask);
scenarioInfo.refScen = refScen;
end

function selectedCtScenIds = resolveScen4DSelection(pln,scenarioCtScenIds,matRad_cfg)
activeCtScenIds = unique(scenarioCtScenIds(:),'stable');
selectedCtScenIds = 1;

if ~isfield(pln,'propOpt') || ~isstruct(pln.propOpt) || ...
        ~isfield(pln.propOpt,'scen4D') || isempty(pln.propOpt.scen4D)
    validateSelectedCtScenIds(selectedCtScenIds,activeCtScenIds,matRad_cfg);
    return;
end

scen4D = pln.propOpt.scen4D;
if isstring(scen4D) && isscalar(scen4D)
    scen4D = char(scen4D);
end

if ischar(scen4D)
    scen4D = strtrim(scen4D);
    if strcmpi(scen4D,'all')
        selectedCtScenIds = activeCtScenIds;
        return;
    end
    matRad_cfg.dispError(['pln.propOpt.scen4D must be ''all'' or a numeric ', ...
        'vector of positive integer CT scenario ids.']);
end

if ~isnumeric(scen4D) || isempty(scen4D) || any(~isfinite(scen4D(:))) || ...
        any(scen4D(:) < 1) || any(fix(scen4D(:)) ~= scen4D(:))
    matRad_cfg.dispError(['pln.propOpt.scen4D must be ''all'' or a numeric ', ...
        'vector of positive integer CT scenario ids.']);
end

selectedCtScenIds = unique(double(scen4D(:)),'stable');
validateSelectedCtScenIds(selectedCtScenIds,activeCtScenIds,matRad_cfg);
end

function validateSelectedCtScenIds(selectedCtScenIds,activeCtScenIds,matRad_cfg)
missingCtScenIds = selectedCtScenIds(~ismember(selectedCtScenIds,activeCtScenIds));
if ~isempty(missingCtScenIds)
    matRad_cfg.dispError('Scenario dose scen4D selection includes inactive CT scenario ids: %s.', ...
        mat2str(missingCtScenIds(:)'));
end
end

function value = validatePositiveInteger(value,name,matRad_cfg)
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value < 1 || fix(value) ~= value
    matRad_cfg.dispError('%s must be a positive integer scalar.',name);
end
value = double(value);
end
