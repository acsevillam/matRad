function scenarioMaps = matRad_buildScenarioDoseMappings(ct, dij, scenarioCtScenIds, refScen, matRadCfg)
% matRad_buildScenarioDoseMappings builds reference CT mapping metadata
%
% call
%   scenarioMaps = ScenarioBatch.Mapping.matRad_buildScenarioDoseMappings(ct,dij,scenarioCtScenIds,refScen,matRadCfg)
%
% input
%   ct:                matRad ct struct; non-reference CT scenarios require
%                      pull DVFs in ct.dvf
%   dij:               matRad dose influence struct with dose grid metadata
%   scenarioCtScenIds: CT scenario id for each active dose scenario
%   refScen:           reference CT scenario id used for mapped rows
%   matRadCfg:        MatRad_Config instance for diagnostics
%
% output
%   scenarioMaps:      cell array with one mapping struct per active
%                      scenario; unmapped reference scenarios are marked by
%                      mapped = false
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

scenarioMaps = cell(numel(scenarioCtScenIds), 1);
for s = 1:numel(scenarioCtScenIds)
    if scenarioCtScenIds(s) == refScen
        scenarioMaps{s}.mapped = false;
    else
        scenarioMaps{s} = matRad_buildDvfMapping(ct, dij, scenarioCtScenIds(s), refScen, matRadCfg);
    end
end
end

function map = matRad_buildDvfMapping(ct, dij, ctScenId, refScen, matRadCfg)
matRad_validateReferenceScenarioMetadata(ct, refScen, matRadCfg);

if ~strcmp(matRad_getDvfType(ct), 'pull')
    matRadCfg.dispError('MultiCT scenario dose calculation requires pull DVFs.');
end

dvf = matRad_getScenarioDvf(ct, ctScenId, matRadCfg);
[dvfX, dvfY, dvfZ] = matRad_getDvfComponents(dvf, matRadCfg);
doseDim = matRad_getDoseGridDimensions(dij, matRadCfg);
sourceDim = [size(dvfX, 1) size(dvfX, 2) size(dvfX, 3)];
dvfUnits = matRad_getDvfUnits(ct);

if isequal(sourceDim, doseDim)
    [dvfX, dvfY, dvfZ] = matRad_convertDvfToDoseVoxelUnits(dvfX, dvfY, dvfZ, ...
                                                           ct, dij, dvfUnits, true, matRadCfg);
elseif isfield(ct, 'cubeDim') && isequal(sourceDim, ct.cubeDim)
    [dvfX, dvfY, dvfZ] = matRad_convertDvfToMillimeters(dvfX, dvfY, dvfZ, ct, dvfUnits, matRadCfg);
    [dvfX, dvfY, dvfZ] = matRad_resampleDvfToDoseGrid(dvfX, dvfY, dvfZ, ct, dij, matRadCfg);
    [dvfX, dvfY, dvfZ] = matRad_millimetersToDoseVoxels(dvfX, dvfY, dvfZ, dij, matRadCfg);
else
    matRadCfg.dispError('DVF dimensions must match ct.cubeDim or dij.doseGrid.dimensions.');
end

map.mapped = true;
map.sourceCtScenId = ctScenId;
map.referenceCtScenId = refScen;
map.dvfX = dvfX;
map.dvfY = dvfY;
map.dvfZ = dvfZ;
map.doseDim = doseDim;
end

function matRad_validateReferenceScenarioMetadata(ct, refScen, matRadCfg)
if isfield(ct, 'refScen') && ~isempty(ct.refScen) && ct.refScen ~= refScen
    matRadCfg.dispError('Requested refScen %d does not match ct.refScen (%d).', ...
                        refScen, ct.refScen);
end

if isfield(ct, 'dvfMetadata')
    metadataRef = [];
    if isfield(ct.dvfMetadata, 'refScen') && ~isempty(ct.dvfMetadata.refScen)
        metadataRef = ct.dvfMetadata.refScen;
    elseif isfield(ct.dvfMetadata, 'referenceCtScen') && ~isempty(ct.dvfMetadata.referenceCtScen)
        metadataRef = ct.dvfMetadata.referenceCtScen;
    end
    if ~isempty(metadataRef) && metadataRef ~= refScen
        matRadCfg.dispError('Requested refScen %d does not match DVF reference scenario %d.', ...
                            refScen, metadataRef);
    end
end
end

function dvfType = matRad_getDvfType(ct)
dvfType = '';
if isfield(ct, 'dvfMetadata') && isfield(ct.dvfMetadata, 'dvfType') && ...
   ~isempty(ct.dvfMetadata.dvfType)
    dvfType = char(ct.dvfMetadata.dvfType);
elseif isfield(ct, 'dvfType') && ~isempty(ct.dvfType)
    dvfType = char(ct.dvfType);
end
end

function dvfUnits = matRad_getDvfUnits(ct)
dvfUnits = 'mm';
if isfield(ct, 'dvfMetadata')
    if isfield(ct.dvfMetadata, 'dvfUnits') && ~isempty(ct.dvfMetadata.dvfUnits)
        dvfUnits = char(ct.dvfMetadata.dvfUnits);
    elseif isfield(ct.dvfMetadata, 'units') && ~isempty(ct.dvfMetadata.units)
        dvfUnits = char(ct.dvfMetadata.units);
    end
end
if strcmp(dvfUnits, 'mm') && isfield(ct, 'dvfUnits') && ~isempty(ct.dvfUnits)
    dvfUnits = char(ct.dvfUnits);
end
end

function dvf = matRad_getScenarioDvf(ct, ctScenId, matRadCfg)
if ~isfield(ct, 'dvf') || isempty(ct.dvf) || numel(ct.dvf) < ctScenId || ...
   isempty(ct.dvf{ctScenId})
    matRadCfg.dispError('ct.dvf must contain a DVF for CT scenario %d.', ctScenId);
end
dvf = ct.dvf{ctScenId};
end

function [dvfX, dvfY, dvfZ] = matRad_getDvfComponents(dvf, matRadCfg)
if size(dvf, 1) == 3
    dvfX = reshape(dvf(1, :, :, :), size(dvf, 2), size(dvf, 3), size(dvf, 4));
    dvfY = reshape(dvf(2, :, :, :), size(dvf, 2), size(dvf, 3), size(dvf, 4));
    dvfZ = reshape(dvf(3, :, :, :), size(dvf, 2), size(dvf, 3), size(dvf, 4));
elseif ndims(dvf) == 4 && size(dvf, 4) == 3
    dvfX = reshape(dvf(:, :, :, 1), size(dvf, 1), size(dvf, 2), size(dvf, 3));
    dvfY = reshape(dvf(:, :, :, 2), size(dvf, 1), size(dvf, 2), size(dvf, 3));
    dvfZ = reshape(dvf(:, :, :, 3), size(dvf, 1), size(dvf, 2), size(dvf, 3));
else
    matRadCfg.dispError('DVF must store three vector components.');
end
end

function doseDim = matRad_getDoseGridDimensions(dij, matRadCfg)
if ~isfield(dij, 'doseGrid') || ~isfield(dij.doseGrid, 'dimensions')
    matRadCfg.dispError('dij.doseGrid.dimensions is required.');
end
doseDim = dij.doseGrid.dimensions;
end

function [dvfX, dvfY, dvfZ] = matRad_convertDvfToDoseVoxelUnits(dvfX, dvfY, dvfZ, ct, dij, dvfUnits, isDoseGrid, matRadCfg)
if strcmpi(dvfUnits, 'voxel') || strcmpi(dvfUnits, 'voxels') || strcmpi(dvfUnits, 'grid')
    if isDoseGrid
        return
    end
end

[dvfX, dvfY, dvfZ] = matRad_convertDvfToMillimeters(dvfX, dvfY, dvfZ, ct, dvfUnits, matRadCfg);
[dvfX, dvfY, dvfZ] = matRad_millimetersToDoseVoxels(dvfX, dvfY, dvfZ, dij, matRadCfg);
end

function [dvfX, dvfY, dvfZ] = matRad_convertDvfToMillimeters(dvfX, dvfY, dvfZ, ct, dvfUnits, matRadCfg)
switch lower(dvfUnits)
    case {'mm', 'millimeter', 'millimeters'}
        % Already in millimeters.
    case {'voxel', 'voxels', 'grid'}
        ctResolution = matRad_getGridResolution(ct, matRadCfg);
        dvfX = dvfX .* ctResolution.x;
        dvfY = dvfY .* ctResolution.y;
        dvfZ = dvfZ .* ctResolution.z;
    otherwise
        matRadCfg.dispError('Unsupported DVF units: %s.', dvfUnits);
end
end

function [dvfX, dvfY, dvfZ] = matRad_millimetersToDoseVoxels(dvfX, dvfY, dvfZ, dij, matRadCfg)
doseResolution = matRad_getGridResolution(dij.doseGrid, matRadCfg);
dvfX = dvfX ./ doseResolution.x;
dvfY = dvfY ./ doseResolution.y;
dvfZ = dvfZ ./ doseResolution.z;
end

function resolution = matRad_getGridResolution(gridStruct, matRadCfg)
if isfield(gridStruct, 'resolution')
    resolution = gridStruct.resolution;
    return
end

if all(isfield(gridStruct, {'x', 'y', 'z'}))
    resolution.x = mean(diff(gridStruct.x));
    resolution.y = mean(diff(gridStruct.y));
    resolution.z = mean(diff(gridStruct.z));
else
    matRadCfg.dispError('Grid resolution or x/y/z vectors are required for DVF unit conversion.');
end
end

function [dvfX, dvfY, dvfZ] = matRad_resampleDvfToDoseGrid(dvfX, dvfY, dvfZ, ct, dij, matRadCfg)
[ctX, ctY, ctZ, doseX, doseY, doseZ] = matRad_getResamplingAxes(ct, dij, ...
                                                                [size(dvfX, 1) size(dvfX, 2) size(dvfX, 3)], matRadCfg);

dvfX = matRad_interp3(ctX, ctY', ctZ, dvfX, doseX, doseY', doseZ, 'linear', 0);
dvfY = matRad_interp3(ctX, ctY', ctZ, dvfY, doseX, doseY', doseZ, 'linear', 0);
dvfZ = matRad_interp3(ctX, ctY', ctZ, dvfZ, doseX, doseY', doseZ, 'linear', 0);
end

function [ctX, ctY, ctZ, doseX, doseY, doseZ] = matRad_getResamplingAxes(ct, dij, sourceDim, matRadCfg)
[ctX, ctY, ctZ] = matRad_getCtAxes(ct, dij, matRadCfg);
[doseX, doseY, doseZ] = matRad_getDoseAxes(dij, matRadCfg);

if numel(ctY) ~= sourceDim(1) || numel(ctX) ~= sourceDim(2) || ...
   numel(ctZ) ~= sourceDim(3)
    matRadCfg.dispError('CT grid axes are inconsistent with DVF dimensions.');
end
end

function [ctX, ctY, ctZ] = matRad_getCtAxes(ct, dij, matRadCfg)
if all(isfield(ct, {'x', 'y', 'z'})) && ~isempty(ct.x) && ~isempty(ct.y) && ~isempty(ct.z)
    ctX = ct.x;
    ctY = ct.y;
    ctZ = ct.z;
    return
end

if isfield(dij, 'ctGrid') && all(isfield(dij.ctGrid, {'x', 'y', 'z'})) && ...
   ~isempty(dij.ctGrid.x) && ~isempty(dij.ctGrid.y) && ~isempty(dij.ctGrid.z)
    ctX = dij.ctGrid.x;
    ctY = dij.ctGrid.y;
    ctZ = dij.ctGrid.z;
    return
end

if isfield(ct, 'cubeDim') && isfield(ct, 'resolution')
    ctWithAxes = matRad_getWorldAxes(ct);
    ctX = ctWithAxes.x;
    ctY = ctWithAxes.y;
    ctZ = ctWithAxes.z;
    return
end

matRadCfg.dispError('CT axes are required to resample DVFs to the dose grid.');
end

function [doseX, doseY, doseZ] = matRad_getDoseAxes(dij, matRadCfg)
if isfield(dij, 'doseGrid') && all(isfield(dij.doseGrid, {'x', 'y', 'z'})) && ...
   ~isempty(dij.doseGrid.x) && ~isempty(dij.doseGrid.y) && ~isempty(dij.doseGrid.z)
    doseX = dij.doseGrid.x;
    doseY = dij.doseGrid.y;
    doseZ = dij.doseGrid.z;
    return
end

if isfield(dij, 'doseGrid') && isfield(dij.doseGrid, 'dimensions') && ...
   isfield(dij.doseGrid, 'resolution')
    doseGrid = matRad_getWorldAxes(dij.doseGrid);
    doseX = doseGrid.x;
    doseY = doseGrid.y;
    doseZ = doseGrid.z;
    return
end

matRadCfg.dispError('Dose grid axes are required to resample DVFs.');
end
