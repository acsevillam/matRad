function [doseCubeRef,meta] = matRad_mapDoseToReferenceScenario(ct,doseCube,ctScenId,refScen)
% matRad_mapDoseToReferenceScenario maps a dose cube to the reference CT scenario
%
% call
%   doseCubeRef = matRad_mapDoseToReferenceScenario(ct,doseCube,ctScenId)
%   [doseCubeRef,meta] = matRad_mapDoseToReferenceScenario(ct,doseCube,ctScenId,refScen)
%
% input
%   ct:          matRad ct struct with pull deformation vector fields
%   doseCube:    dose cube calculated on CT scenario ctScenId
%   ctScenId:    CT scenario id of doseCube
%   refScen:     reference CT scenario id (default: ct.refScen or 1)
%
% output
%   doseCubeRef: doseCube mapped to the reference CT scenario
%   meta:        mapping metadata
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

matRad_cfg = MatRad_Config.instance();

if nargin < 4 || isempty(refScen)
    refScen = getReferenceScenario(ct);
end

ctScenId = validateCtScenarioId(ctScenId,'ctScenId',ct,matRad_cfg);
refScen = validateCtScenarioId(refScen,'refScen',ct,matRad_cfg);
validateReferenceScenarioMetadata(ct,refScen,matRad_cfg);

if ~isequal(size(doseCube),ct.cubeDim)
    matRad_cfg.dispError('doseCube size must match ct.cubeDim.');
end

meta = struct();
meta.sourceCtScenId = ctScenId;
meta.referenceCtScenId = refScen;
meta.dvfType = getDvfType(ct);
meta.dvfUnits = getDvfUnits(ct);
meta.method = 'pullDirectDoseMapping';
meta.mapped = ctScenId ~= refScen;

if ctScenId == refScen
    doseCubeRef = doseCube;
    return;
end

if ~strcmp(meta.dvfType,'pull')
    matRad_cfg.dispError('Mapping dose to a reference CT scenario requires pull DVFs.');
end

dvf = getScenarioDvf(ct,ctScenId,matRad_cfg);
[dvfX,dvfY,dvfZ] = getDvfComponents(dvf,ct,matRad_cfg);
[dvfX,dvfY,dvfZ] = convertDvfToVoxelUnits(dvfX,dvfY,dvfZ,ct,meta.dvfUnits,matRad_cfg);

xGridVec = 1:ct.cubeDim(2);
yGridVec = 1:ct.cubeDim(1);
zGridVec = 1:ct.cubeDim(3);
[X,Y,Z] = meshgrid(xGridVec,yGridVec,zGridVec);

doseCubeRef = matRad_interp3(X,Y,Z,doseCube, ...
    X - dvfX,Y - dvfY,Z - dvfZ,'linear',0);

end

function refScen = getReferenceScenario(ct)
if isfield(ct,'refScen') && ~isempty(ct.refScen)
    refScen = ct.refScen;
else
    refScen = 1;
end
end

function scen = validateCtScenarioId(scen,name,ct,matRad_cfg)
if ~(isnumeric(scen) && isscalar(scen) && isfinite(scen) && round(scen) == scen && scen >= 1)
    matRad_cfg.dispError('%s must be a positive integer scalar.',name);
end
scen = double(scen);

if isfield(ct,'numOfCtScen') && scen > ct.numOfCtScen
    matRad_cfg.dispError('%s (%d) exceeds ct.numOfCtScen (%d).',name,scen,ct.numOfCtScen);
end
end

function validateReferenceScenarioMetadata(ct,refScen,matRad_cfg)
ctRefScen = getOptionalReferenceScenario(ct);
if ~isempty(ctRefScen) && ctRefScen ~= refScen
    matRad_cfg.dispError('Requested reference CT scenario %d does not match ct.refScen (%d).', ...
        refScen,ctRefScen);
end

dvfRefScen = getOptionalDvfReferenceScenario(ct);
if ~isempty(dvfRefScen) && dvfRefScen ~= refScen
    matRad_cfg.dispError('Requested reference CT scenario %d does not match the DVF reference CT scenario (%d).', ...
        refScen,dvfRefScen);
end
end

function refScen = getOptionalReferenceScenario(ct)
refScen = [];
if isfield(ct,'refScen') && ~isempty(ct.refScen)
    refScen = ct.refScen;
end
end

function refScen = getOptionalDvfReferenceScenario(ct)
refScen = [];
if ~isfield(ct,'dvfMetadata')
    return;
end

if isfield(ct.dvfMetadata,'refScen') && ~isempty(ct.dvfMetadata.refScen)
    refScen = ct.dvfMetadata.refScen;
elseif isfield(ct.dvfMetadata,'referenceCtScen') && ~isempty(ct.dvfMetadata.referenceCtScen)
    refScen = ct.dvfMetadata.referenceCtScen;
elseif isfield(ct.dvfMetadata,'referenceScenario') && ~isempty(ct.dvfMetadata.referenceScenario)
    refScen = ct.dvfMetadata.referenceScenario;
end
end

function dvfType = getDvfType(ct)
if isfield(ct,'dvfMetadata') && isfield(ct.dvfMetadata,'dvfType') && ~isempty(ct.dvfMetadata.dvfType)
    dvfType = char(ct.dvfMetadata.dvfType);
elseif isfield(ct,'dvfType') && ~isempty(ct.dvfType)
    dvfType = char(ct.dvfType);
else
    dvfType = '';
end
end

function dvfUnits = getDvfUnits(ct)
dvfUnits = 'mm';
if isfield(ct,'dvfMetadata')
    if isfield(ct.dvfMetadata,'dvfUnits') && ~isempty(ct.dvfMetadata.dvfUnits)
        dvfUnits = char(ct.dvfMetadata.dvfUnits);
    elseif isfield(ct.dvfMetadata,'units') && ~isempty(ct.dvfMetadata.units)
        dvfUnits = char(ct.dvfMetadata.units);
    end
end

if strcmp(dvfUnits,'mm') && isfield(ct,'dvfUnits') && ~isempty(ct.dvfUnits)
    dvfUnits = char(ct.dvfUnits);
end
end

function dvf = getScenarioDvf(ct,ctScenId,matRad_cfg)
if ~isfield(ct,'dvf') || isempty(ct.dvf) || numel(ct.dvf) < ctScenId || isempty(ct.dvf{ctScenId})
    matRad_cfg.dispError('ct.dvf must contain a deformation vector field for CT scenario %d.',ctScenId);
end
dvf = ct.dvf{ctScenId};
end

function [dvfX,dvfY,dvfZ] = getDvfComponents(dvf,ct,matRad_cfg)
if ndims(dvf) ~= 4
    matRad_cfg.dispError('DVF must be a 4-D array with three vector components.');
end

if size(dvf,1) == 3
    dvfX = squeeze(dvf(1,:,:,:));
    dvfY = squeeze(dvf(2,:,:,:));
    dvfZ = squeeze(dvf(3,:,:,:));
elseif size(dvf,4) == 3
    dvfX = squeeze(dvf(:,:,:,1));
    dvfY = squeeze(dvf(:,:,:,2));
    dvfZ = squeeze(dvf(:,:,:,3));
else
    matRad_cfg.dispError('DVF must store three vector components in the first or fourth dimension.');
end

if ~isequal(size(dvfX),ct.cubeDim) || ~isequal(size(dvfY),ct.cubeDim) || ~isequal(size(dvfZ),ct.cubeDim)
    matRad_cfg.dispError('DVF component size must match ct.cubeDim.');
end
end

function [dvfX,dvfY,dvfZ] = convertDvfToVoxelUnits(dvfX,dvfY,dvfZ,ct,dvfUnits,matRad_cfg)
switch lower(dvfUnits)
    case {'mm','millimeter','millimeters'}
        dvfX = dvfX ./ ct.resolution.x;
        dvfY = dvfY ./ ct.resolution.y;
        dvfZ = dvfZ ./ ct.resolution.z;
    case {'voxel','voxels','grid'}
        % Already in intrinsic voxel coordinates.
    otherwise
        matRad_cfg.dispError('Unsupported DVF units: %s.',dvfUnits);
end
end
