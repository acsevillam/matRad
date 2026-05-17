function [ctOut, commonGrid, cropIdx] = matRad_cropCtToOverlap(ctIn, tol)
% matRad function to crop multiple CT scenarios to their common overlap
%
% call
%   [ctOut,commonGrid,cropIdx] = matRad_cropCtToOverlap(ctIn)
%   [ctOut,commonGrid,cropIdx] = matRad_cropCtToOverlap(ctIn,tol)
%
% input
%   ctIn:           cell array with one matRad ct struct per CT scenario
%   tol:            optional numeric tolerance for axis spacing and overlap
%                   checks in mm
%
% output
%   ctOut:          cell array with CT scenarios cropped to their common
%                   physical overlap
%   commonGrid:     struct with x, y, and z axes of the common cropped grid
%   cropIdx:        cell array with Ix, Iy, and Iz index vectors used to
%                   crop each original CT scenario
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

matRadCfg = MatRad_Config.instance();

if nargin < 2 || isempty(tol)
    tol = 1e-5;
end

if ~iscell(ctIn) || isempty(ctIn)
    matRadCfg.dispError('ctIn must be a non-empty cell array of ct structs.');
end

numOfScenarios = numel(ctIn);

allX = cellfun(@(ct) ct.x, ctIn, 'UniformOutput', false);
allY = cellfun(@(ct) ct.y, ctIn, 'UniformOutput', false);
allZ = cellfun(@(ct) ct.z, ctIn, 'UniformOutput', false);

dx = cellfun(@(axis) median(diff(axis)), allX);
dy = cellfun(@(axis) median(diff(axis)), allY);
dz = cellfun(@(axis) median(diff(axis)), allZ);

if max(dx) - min(dx) > tol
    matRadCfg.dispError('CT scenarios have different x-axis spacing. Resampling is required.');
end
if max(dy) - min(dy) > tol
    matRadCfg.dispError('CT scenarios have different y-axis spacing. Resampling is required.');
end
if max(dz) - min(dz) > tol
    matRadCfg.dispError('CT scenarios have different z-axis spacing. Resampling is required.');
end

xMinCommon = max(cellfun(@(axis) min(axis), allX));
xMaxCommon = min(cellfun(@(axis) max(axis), allX));
yMinCommon = max(cellfun(@(axis) min(axis), allY));
yMaxCommon = min(cellfun(@(axis) max(axis), allY));
zMinCommon = max(cellfun(@(axis) min(axis), allZ));
zMaxCommon = min(cellfun(@(axis) max(axis), allZ));

if ~(xMinCommon < xMaxCommon && yMinCommon < yMaxCommon && zMinCommon < zMaxCommon)
    matRadCfg.dispError('CT scenarios do not share a common 3-D overlap.');
end

ctOut = ctIn;
cropIdx = cell(1, numOfScenarios);

for ctScenario = 1:numOfScenarios
    xAxis = ctIn{ctScenario}.x;
    yAxis = ctIn{ctScenario}.y;
    zAxis = ctIn{ctScenario}.z;

    ixCrop = find(xAxis >= xMinCommon - tol & xAxis <= xMaxCommon + tol);
    iyCrop = find(yAxis >= yMinCommon - tol & yAxis <= yMaxCommon + tol);
    izCrop = find(zAxis >= zMinCommon - tol & zAxis <= zMaxCommon + tol);

    if isempty(ixCrop) || isempty(iyCrop) || isempty(izCrop)
        matRadCfg.dispError('CT scenario %d has no voxels inside the common overlap.', ctScenario);
    end

    ctCube = ctIn{ctScenario}.cubeHU{1};
    croppedCube = ctCube(iyCrop, ixCrop, izCrop);
    croppedDim = matRad_getCubeDim(croppedCube);
    originalDim = matRad_getCubeDim(ctCube);

    ctOut{ctScenario}.x = xAxis(ixCrop);
    ctOut{ctScenario}.y = yAxis(iyCrop);
    ctOut{ctScenario}.z = zAxis(izCrop);
    ctOut{ctScenario}.cubeHU{1} = croppedCube;
    ctOut{ctScenario}.cubeDim = croppedDim;

    cropIdx{ctScenario} = struct('Ix', ixCrop, 'Iy', iyCrop, 'Iz', izCrop, ...
                                 'origDim', originalDim);
end

nx = cellfun(@(ct) numel(ct.x), ctOut);
ny = cellfun(@(ct) numel(ct.y), ctOut);
nz = cellfun(@(ct) numel(ct.z), ctOut);
if ~(all(nx == nx(1)) && all(ny == ny(1)) && all(nz == nz(1)))
    matRadCfg.dispError('CT scenarios still have different dimensions after overlap cropping.');
end

commonGrid = struct('x', ctOut{1}.x, 'y', ctOut{1}.y, 'z', ctOut{1}.z);
end

function cubeDim = matRad_getCubeDim(cube)
cubeDim = size(cube);
if numel(cubeDim) < 3
    cubeDim(3) = 1;
end
cubeDim = double(cubeDim(1:3));
end
