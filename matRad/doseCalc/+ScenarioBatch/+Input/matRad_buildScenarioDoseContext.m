function ctx = matRad_buildScenarioDoseContext(ct, cst, pln, dij, cfg, calculationMode, matRadCfg, scenarioInfo)
% matRad_buildScenarioDoseContext validates shared scenario-dose inputs
%
% call
%   ctx = ScenarioBatch.Input.matRad_buildScenarioDoseContext(ct,cst,pln,dij,cfg,calculationMode,matRadCfg)
%   ctx = ScenarioBatch.Input.matRad_buildScenarioDoseContext(ct,cst,pln,dij,cfg,calculationMode,matRadCfg,scenarioInfo)
%
% input
%   ct:           matRad ct struct
%   cst:          matRad cst cell array
%   pln:          matRad pln struct with pln.multScen as matRad_ScenarioModel
%   dij:          robust dose influence struct containing the requested
%                 linear quantity as scenario cell matrices
%   cfg:          robust scenario-dose configuration struct; dose quantities
%                 are in Gy or Gy(RBE) according to the selected linear dij
%                 field
%   calculationMode: method identifier, e.g. 'PROB', 'INTERVAL2', or
%                    'INTERVAL3'
%   matRadCfg:   MatRad_Config instance for diagnostics
%   scenarioInfo: optional output of
%                 ScenarioBatch.Input.matRad_selectScenarioDoseRealizations
%
% output
%   ctx:          validated context struct with quantity metadata, DIJ
%                 scenario indices, CT scenario ids selected by
%                 pln.propOpt.scen4D (default: 1), normalized scenario
%                 weights, selected target/OAR rows, selected structure
%                 voxels in the reference scenario, and CT mapping metadata
%
% References
%   -
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

cfg = ScenarioBatch.Config.matRad_normalizeScenarioDoseInputConfig(cfg, ct, pln, ...
                                                                   calculationMode, matRadCfg);
quantity = ScenarioBatch.Quantity.matRad_resolveScenarioDoseQuantity(dij, pln, cfg, matRadCfg);
if nargin < 8 || isempty(scenarioInfo)
    scenarioInfo = ScenarioBatch.Input.matRad_selectScenarioDoseRealizations(ct, pln, cfg, matRadCfg);
end
matRad_validateActiveScenarioCells(dij, quantity.field, scenarioInfo.scenarioDijIx, matRadCfg);
scenarioDijIx = scenarioInfo.scenarioDijIx(:);
scenarioCtScenIds = scenarioInfo.scenarioCtScenIds(:);
scenarioWeights = scenarioInfo.scenarioWeights(:);

structureInfo = ScenarioBatch.Structures.matRad_selectScenarioDoseStructures(cst, dij, cfg);

ctx.cfg = cfg;
ctx.quantity = quantity;
ctx.scenarioDijIx = scenarioDijIx;
ctx.scenarioCtScenIds = scenarioCtScenIds;
ctx.scenarioWeights = scenarioWeights;
ctx.targetRows = structureInfo.targetRows;
ctx.oarRows = structureInfo.oarRows;
ctx.targetStructures = structureInfo.targetStructures;
ctx.oarStructures = structureInfo.oarStructures;
ctx.targetStructRows = structureInfo.targetStructRows;
ctx.oarStructRows = structureInfo.oarStructRows;
ctx.cstDoseGrid = structureInfo.cstDoseGrid;
ctx.numVoxels = matRad_getNumDoseVoxels(dij, quantity.matrixCell, matRadCfg);
ctx.numBixels = matRad_getNumBixels(dij, quantity.matrixCell, scenarioDijIx, matRadCfg);
ctx.scenarioMaps = ScenarioBatch.Mapping.matRad_buildScenarioDoseMappings(ct, dij, ...
                                                                          scenarioCtScenIds, cfg.refScen, matRadCfg);
end

function matRad_validateActiveScenarioCells(dij, quantityField, scenarioDijIx, matRadCfg)
if isempty(scenarioDijIx)
    matRadCfg.dispError('No active scenarios found for scenario dose calculation.');
end

quantityCells = dij.(quantityField);
if max(scenarioDijIx) > numel(quantityCells)
    matRadCfg.dispError('Scenario indices exceed dij.%s cell array dimensions.', quantityField);
end

emptyScenario = cellfun(@isempty, quantityCells(scenarioDijIx));
if any(emptyScenario)
    matRadCfg.dispError('dij.%s contains empty active scenarios.', quantityField);
end
end

function numVoxels = matRad_getNumDoseVoxels(dij, matrixCell, matRadCfg)
if isfield(dij, 'doseGrid') && isfield(dij.doseGrid, 'numOfVoxels')
    numVoxels = dij.doseGrid.numOfVoxels;
elseif isfield(dij, 'doseGrid') && isfield(dij.doseGrid, 'dimensions')
    numVoxels = prod(dij.doseGrid.dimensions);
else
    firstNonEmpty = find(~cellfun(@isempty, matrixCell(:)), 1, 'first');
    if isempty(firstNonEmpty)
        matRadCfg.dispError('Could not determine number of dose voxels.');
    end
    numVoxels = size(matrixCell{firstNonEmpty}, 1);
end
end

function numBixels = matRad_getNumBixels(dij, matrixCell, scenarioDijIx, matRadCfg)
if isfield(dij, 'totalNumOfBixels') && ~isempty(dij.totalNumOfBixels)
    numBixels = dij.totalNumOfBixels;
else
    numBixels = size(matrixCell{scenarioDijIx(1)}, 2);
end

if any(cellfun(@(m) size(m, 2) ~= numBixels, matrixCell(scenarioDijIx)))
    matRadCfg.dispError('All active scenario matrices must have the same number of bixels.');
end
end
