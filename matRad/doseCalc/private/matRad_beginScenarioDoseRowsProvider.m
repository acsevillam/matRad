function [provider,source] = matRad_beginScenarioDoseRowsProvider(provider,ctx,quantity,scenarioRowIx,matRad_cfg)
% matRad_beginScenarioDoseRowsProvider prepares one scenario row source
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
if strcmp(provider.type,'inMemory')
    source.quantity = quantity;
    source.scenarioDijIx = ctx.scenarioDijIx(scenarioRowIx);
    return;
end

scenarioId = provider.scenarioIds(scenarioRowIx);
if isfield(provider,'preloadedScenarioId') && ...
        isequal(provider.preloadedScenarioId,scenarioId)
    dijScenario = provider.preloadedDij;
    provider = rmfield(provider,{'preloadedScenarioId','preloadedDij'});
else
    pln_s = matRad_prepareSerialScenarioPlan(provider.pln,scenarioId);
    dijScenario = matRad_calcDoseInfluence(provider.ct,provider.cst, ...
        provider.stf,pln_s);
end

matrix = extractSingleScenarioQuantityMatrix(dijScenario,quantity.field, ...
    ctx.numBixels,matRad_cfg);
singleQuantity = quantity;
singleQuantity.matrixCell = {matrix};
source.quantity = singleQuantity;
source.scenarioDijIx = 1;
end

function matrix = extractSingleScenarioQuantityMatrix(dijScenario,fieldName,numBixels,matRad_cfg)
if ~isfield(dijScenario,fieldName) || ~iscell(dijScenario.(fieldName))
    matRad_cfg.dispError('Single-scenario dij.%s is missing or is not a cell array.',fieldName);
end

quantityCells = dijScenario.(fieldName);
activeIx = find(~cellfun(@isempty,quantityCells(:)));
if isempty(activeIx)
    matRad_cfg.dispError('Single-scenario dij.%s contains no active scenario matrix.',fieldName);
end

matrix = quantityCells{activeIx(1)};
if size(matrix,2) ~= numBixels
    matRad_cfg.dispError('Single-scenario dij.%s bixel count is inconsistent with the streaming context.',fieldName);
end
end
