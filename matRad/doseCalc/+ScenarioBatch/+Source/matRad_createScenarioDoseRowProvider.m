function [provider, dijForResolve] = matRad_createScenarioDoseRowProvider( ...
                                                                          ct, cst, stf, pln, cfg, scenarioInfo, matRadCfg, calculationName)
% matRad_createScenarioDoseRowProvider creates scenario-batch row access
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

if nargin < 8 || isempty(calculationName)
    calculationName = 'scenario-batch scenario dose';
end

provider = struct();
provider.secondPassStrategy = cfg.SecondPassStrategy;
provider.keepCache = cfg.KeepCache;
provider.cacheRoot = cfg.CacheRoot;
provider.cacheContext = [];

if isfield(cfg, 'PrecomputedDij') && ~isempty(cfg.PrecomputedDij)
    provider.type = 'inMemory';
    provider.dij = cfg.PrecomputedDij;
    dijForResolve = cfg.PrecomputedDij;
    return
end

provider.type = 'scenario-batch';
provider.ct = ct;
provider.cst = cst;
provider.stf = stf;
provider.pln = ScenarioBatch.Worker.matRad_sanitizeWorkerPlan(pln);
provider.scenarioIds = scenarioInfo.scenarioIds(:);

firstScenarioId = scenarioInfo.scenarioIds(1);
pln_s = ScenarioBatch.Scenarios.matRad_buildSingleScenarioPlan(provider.pln, firstScenarioId);
dijTemplate = matRad_calcDoseInfluence(ct, cst, stf, pln_s);
provider.preloadedScenarioId = firstScenarioId;
provider.preloadedDij = dijTemplate;

dijForResolve = matRad_buildSyntheticDijForScenarioResolution( ...
                                                              dijTemplate, pln, scenarioInfo.scenarioDijIx, matRadCfg, calculationName);
end

function dij = matRad_buildSyntheticDijForScenarioResolution(templateDij, pln, scenarioDijIx, matRadCfg, calculationName)
if nargin < 5 || isempty(calculationName)
    calculationName = 'scenario-batch scenario dose';
end

dij = templateDij;
if ~isfield(templateDij, 'doseGrid') || ~isfield(templateDij.doseGrid, 'numOfVoxels') || ...
        ~isfield(templateDij, 'totalNumOfBixels')
    matRadCfg.dispError('%s template dij is missing doseGrid or totalNumOfBixels metadata.', ...
                        calculationName);
end

containerSize = pln.multScen.getDijContainerSize();
fields = fieldnames(templateDij);
for f = 1:numel(fields)
    fieldName = fields{f};
    if ~iscell(templateDij.(fieldName))
        continue
    end

    sample = matRad_getFirstNonEmptyCell(templateDij.(fieldName));
    if isempty(sample) || ~isnumeric(sample) || ~ismatrix(sample)
        continue
    end

    cells = cell(containerSize);
    for ix = scenarioDijIx(:)'
        cells{ix} = spalloc(templateDij.doseGrid.numOfVoxels, ...
                            templateDij.totalNumOfBixels, 0);
    end
    dij.(fieldName) = cells;
end
end

function value = matRad_getFirstNonEmptyCell(cells)
value = [];
firstIx = find(~cellfun(@isempty, cells(:)), 1, 'first');
if ~isempty(firstIx)
    value = cells{firstIx};
end
end
