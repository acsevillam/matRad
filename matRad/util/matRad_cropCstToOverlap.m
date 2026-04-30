function cstOut = matRad_cropCstToOverlap(cstIn,cropIdx,ctOriginal,ctCropped)
% matRad function to crop structure masks to the common CT overlap
%
% call
%   cstOut = matRad_cropCstToOverlap(cstIn,cropIdx,ctOriginal,ctCropped)
%
% input
%   cstIn:          cell array with one matRad cst per CT scenario. The
%                   structure indices are expected in cstIn{ctScenario}{:,4}{1}
%   cropIdx:        cell array with Ix, Iy, and Iz index vectors returned by
%                   matRad_cropCtToOverlap
%   ctOriginal:     cell array with the original uncropped ct structs
%   ctCropped:      cell array with the cropped ct structs
%
% output
%   cstOut:         cell array with one cropped matRad cst per CT scenario
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

numOfScenarios = numel(ctOriginal);
if numOfScenarios == 0
    matRad_cfg.dispError('ctOriginal must contain at least one CT scenario.');
end

if numel(ctCropped) ~= numOfScenarios || numel(cropIdx) ~= numOfScenarios || ...
        numel(cstIn) ~= numOfScenarios
    matRad_cfg.dispError('cstIn, cropIdx, ctOriginal, and ctCropped must contain the same number of scenarios.');
end

cstOut = cstIn;

for ctScenario = 1:numOfScenarios
    cstScenario = cstOut{ctScenario};
    if ~iscell(cstScenario)
        matRad_cfg.dispError('cstIn{%d} must be a matRad cst cell array.',ctScenario);
    end

    originalDim = validateCtDimensions(ctOriginal{ctScenario},'ctOriginal',ctScenario,matRad_cfg);
    croppedDim = validateCtDimensions(ctCropped{ctScenario},'ctCropped',ctScenario,matRad_cfg);
    [iyCrop,ixCrop,izCrop] = validateCropIndices(cropIdx{ctScenario},originalDim,ctScenario,matRad_cfg);

    mapX = zeros(1,originalDim(2),'uint32');
    mapY = zeros(1,originalDim(1),'uint32');
    mapZ = zeros(1,originalDim(3),'uint32');
    mapX(ixCrop) = uint32(1:numel(ixCrop));
    mapY(iyCrop) = uint32(1:numel(iyCrop));
    mapZ(izCrop) = uint32(1:numel(izCrop));

    for structureIx = 1:size(cstScenario,1)
        if size(cstScenario,2) < 4 || ~iscell(cstScenario{structureIx,4}) || ...
                isempty(cstScenario{structureIx,4})
            matRad_cfg.dispError('cstIn{%d}{%d,4} must contain a cell array of structure indices.', ...
                ctScenario,structureIx);
        end

        structureIndices = cstScenario{structureIx,4}{1};
        if isempty(structureIndices)
            cstScenario{structureIx,4}{1} = [];
            continue;
        end

        structureIndices = normalizeStructureIndices(structureIndices,originalDim, ...
            ctScenario,structureIx,matRad_cfg);
        [iy,ix,iz] = ind2sub(originalDim,structureIndices);
        inCrop = mapY(iy) > 0 & mapX(ix) > 0 & mapZ(iz) > 0;

        if ~any(inCrop)
            cstScenario{structureIx,4}{1} = [];
            continue;
        end

        iyCropped = double(mapY(iy(inCrop)));
        ixCropped = double(mapX(ix(inCrop)));
        izCropped = double(mapZ(iz(inCrop)));
        cstScenario{structureIx,4}{1} = sub2ind(croppedDim,iyCropped,ixCropped,izCropped);
        cstScenario{structureIx,4}{1} = cstScenario{structureIx,4}{1}(:);
    end

    cstOut{ctScenario} = cstScenario;
end

end

function cubeDim = validateCtDimensions(ct,ctName,ctScenario,matRad_cfg)
if ~isstruct(ct) || ~isfield(ct,'cubeHU') || ~iscell(ct.cubeHU) || ...
        isempty(ct.cubeHU) || isempty(ct.cubeHU{1})
    matRad_cfg.dispError('%s{%d} must contain cubeHU{1}.',ctName,ctScenario);
end

cubeDim = size(ct.cubeHU{1});
if numel(cubeDim) < 3
    cubeDim(3) = 1;
elseif numel(cubeDim) > 3
    matRad_cfg.dispError('%s{%d}.cubeHU{1} must be a 3-D cube.',ctName,ctScenario);
end

if ~isfield(ct,'cubeDim') || ~isnumeric(ct.cubeDim) || numel(ct.cubeDim) ~= 3
    matRad_cfg.dispError('%s{%d}.cubeDim must be a numeric [numY numX numZ] vector.',ctName,ctScenario);
end

ctCubeDim = double(ct.cubeDim(:)');
cubeDim = double(cubeDim);
if ~isequal(ctCubeDim,cubeDim)
    matRad_cfg.dispError('%s{%d}.cubeDim does not match size(%s{%d}.cubeHU{1}).', ...
        ctName,ctScenario,ctName,ctScenario);
end
end

function [iyCrop,ixCrop,izCrop] = validateCropIndices(cropInfo,originalDim,ctScenario,matRad_cfg)
requiredFields = {'Ix','Iy','Iz'};
for fieldIx = 1:numel(requiredFields)
    if ~isstruct(cropInfo) || ~isfield(cropInfo,requiredFields{fieldIx})
        matRad_cfg.dispError('cropIdx{%d} must contain Ix, Iy, and Iz fields.',ctScenario);
    end
end

ixCrop = validateIndexVector(cropInfo.Ix,originalDim(2),'Ix',ctScenario,matRad_cfg);
iyCrop = validateIndexVector(cropInfo.Iy,originalDim(1),'Iy',ctScenario,matRad_cfg);
izCrop = validateIndexVector(cropInfo.Iz,originalDim(3),'Iz',ctScenario,matRad_cfg);
end

function indexVector = validateIndexVector(indexVector,maxIndex,fieldName,ctScenario,matRad_cfg)
if isempty(indexVector) || ~isnumeric(indexVector) || ~isvector(indexVector) || ...
        any(indexVector(:) < 1) || any(indexVector(:) > maxIndex) || ...
        any(indexVector(:) ~= fix(indexVector(:)))
    matRad_cfg.dispError('cropIdx{%d}.%s must be an integer vector within the original CT dimensions.', ...
        ctScenario,fieldName);
end

indexVector = indexVector(:)';
end

function structureIndices = normalizeStructureIndices(structureIndices,originalDim,ctScenario,structureIx,matRad_cfg)
if islogical(structureIndices)
    if ~isequal(size(structureIndices),originalDim)
        matRad_cfg.dispError('Logical mask for cstIn{%d}{%d,4}{1} does not match the original CT dimensions.', ...
            ctScenario,structureIx);
    end
    structureIndices = find(structureIndices);
elseif isnumeric(structureIndices) && isvector(structureIndices)
    structureIndices = structureIndices(:);
else
    matRad_cfg.dispError('cstIn{%d}{%d,4}{1} must be a numeric index vector or logical mask.', ...
        ctScenario,structureIx);
end

if any(structureIndices < 1) || any(structureIndices > prod(originalDim)) || ...
        any(structureIndices ~= fix(structureIndices))
    matRad_cfg.dispError('cstIn{%d}{%d,4}{1} contains indices outside the original CT dimensions.', ...
        ctScenario,structureIx);
end
end
