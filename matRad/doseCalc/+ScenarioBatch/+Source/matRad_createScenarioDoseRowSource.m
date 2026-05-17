function [provider, source] = matRad_createScenarioDoseRowSource(provider, ctx, quantity, scenarioRowIx, matRadCfg)
% matRad_createScenarioDoseRowSource creates row access for one scenario
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

source = struct();
if strcmp(provider.type, 'inMemory')
    source.quantity = quantity;
    source.scenarioDijIx = ctx.scenarioDijIx(scenarioRowIx);
    return
end

scenarioId = provider.scenarioIds(scenarioRowIx);
if isfield(provider, 'preloadedScenarioId') && ...
        isequal(provider.preloadedScenarioId, scenarioId)
    dijScenario = provider.preloadedDij;
    provider = rmfield(provider, {'preloadedScenarioId', 'preloadedDij'});
else
    pln_s = ScenarioBatch.Scenarios.matRad_buildSingleScenarioPlan(provider.pln, scenarioId);
    dijScenario = matRad_calcDoseInfluence(provider.ct, provider.cst, ...
                                           provider.stf, pln_s);
end

matrix = matRad_extractSingleScenarioQuantityMatrix(dijScenario, quantity.field, ...
                                                    ctx.numBixels, matRadCfg);
singleQuantity = quantity;
singleQuantity.matrixCell = {matrix};
source.quantity = singleQuantity;
source.scenarioDijIx = 1;
end

function matrix = matRad_extractSingleScenarioQuantityMatrix(dijScenario, fieldName, numBixels, matRadCfg)
if ~isfield(dijScenario, fieldName) || ~iscell(dijScenario.(fieldName))
    matRadCfg.dispError('Single-scenario dij.%s is missing or is not a cell array.', fieldName);
end

quantityCells = dijScenario.(fieldName);
activeIx = find(~cellfun(@isempty, quantityCells(:)));
if isempty(activeIx)
    matRadCfg.dispError('Single-scenario dij.%s contains no active scenario matrix.', fieldName);
end

matrix = quantityCells{activeIx(1)};
if size(matrix, 2) ~= numBixels
    matRadCfg.dispError('Single-scenario dij.%s bixel count is inconsistent with the scenario-batch context.', fieldName);
end
end
