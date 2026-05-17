function structureInfo = matRad_selectScenarioDoseStructures(cst, dij, cfg)
% matRad_selectScenarioDoseStructures selects scenario-dose structure rows
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

cstDoseGrid = matRad_normalizeCstForScenarioDoseRows(cst, dij);
[targetRows, targetStructures] = matRad_selectStructureRows(cstDoseGrid, ...
                                                            cfg.targetStructSel, 'TARGET', cfg.refScen);
[oarRows, oarStructures] = matRad_selectStructureRows(cstDoseGrid, ...
                                                      cfg.OARStructSel, 'OAR', cfg.refScen);

structureInfo.cstDoseGrid = cstDoseGrid;
structureInfo.targetRows = targetRows;
structureInfo.oarRows = oarRows;
structureInfo.targetStructures = targetStructures(:);
structureInfo.oarStructures = oarStructures(:);
structureInfo.targetStructRows = matRad_getSelectedStructureRows(targetStructures);
structureInfo.oarStructRows = matRad_getSelectedStructureRows(oarStructures);
end

function cst = matRad_normalizeCstForScenarioDoseRows(cst, dij)
cst = matRad_applyOverlapPrioritiesIfAvailable(cst);
cst = matRad_resizeCstToDoseGrid(cst, dij);
end

function cst = matRad_applyOverlapPrioritiesIfAvailable(cst)
hasPriority = cellfun(@(voiProperties) isstruct(voiProperties) && ...
                      isfield(voiProperties, 'Priority'), cst(:, 5));
if all(hasPriority)
    cst = matRad_setOverlapPriorities(cst);
end
end

function cst = matRad_resizeCstToDoseGrid(cst, dij)
requiredGridFields = {'x', 'y', 'z'};
if isfield(dij, 'ctGrid') && isfield(dij, 'doseGrid') && ...
        all(isfield(dij.ctGrid, requiredGridFields)) && ...
        all(isfield(dij.doseGrid, requiredGridFields))
    cst = matRad_resizeCstToGrid(cst, dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, ...
                                 dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);
end
end

function [rows, structures] = matRad_selectStructureRows(cst, structureSelection, structureType, refScen)
if isempty(structureSelection)
    cstRows = matRad_resolveStructureSelection(cst, 'all', [], structureType);
else
    cstRows = matRad_resolveStructureSelection(cst, 'include', structureSelection, structureType);
end

rows = [];
structures = struct('cstRow', cell(numel(cstRows), 1), ...
                    'type', cell(numel(cstRows), 1), 'voxelIx', cell(numel(cstRows), 1));
for k = 1:numel(cstRows)
    rowIx = cstRows(k);
    voxels = unique(matRad_getStructureVoxelIndices(cst, rowIx, refScen), 'stable');
    structures(k).cstRow = rowIx;
    structures(k).type = structureType;
    structures(k).voxelIx = voxels(:);
    rows = [rows; voxels(:)]; %#ok<AGROW>
end
rows = unique(rows(:), 'stable');
end

function cstRows = matRad_getSelectedStructureRows(structures)
if isempty(structures)
    cstRows = zeros(0, 1);
else
    cstRows = [structures.cstRow]';
end
end

function voxels = matRad_getStructureVoxelIndices(cst, rowIx, refScen)
voxels = [];
if size(cst, 2) < 4 || isempty(cst{rowIx, 4})
    return
end

contours = cst{rowIx, 4};
if iscell(contours) && numel(contours) >= refScen && ~isempty(contours{refScen})
    voxels = contours{refScen}(:);
elseif iscell(contours) && ~isempty(contours{1})
    voxels = contours{1}(:);
elseif isnumeric(contours)
    voxels = contours(:);
end
end
