function [gammaCube,gammaPassRateCell,gammaFig] = matRad_gammaIndex(cube1,cube2,resolution,criteria,slice,n,localglobal,cst,ct)
% gamma index calculation
% according to http://www.ncbi.nlm.nih.gov/pubmed/9608475
%
% call
%   [gammaCube,gammaPassRateCell] = matRad_gammaIndex(cube1,cube2,resolution,criteria,cst)
%   [gammaCube,gammaPassRateCell] = matRad_gammaIndex(cube1,cube2,resolution,criteria,slice,cst)
%   [gammaCube,gammaPassRateCell] = matRad_gammaIndex(cube1,cube2,resolution,criteria,n,cst)
%   [gammaCube,gammaPassRateCell] = matRad_gammaIndex(cube1,cube2,resolution,criteria,localglobal,cst)
%   [gammaCube,gammaPassRateCell] = matRad_gammaIndex(cube1,cube2,resolution,criteria,slice,n,cst)
%               ...
%   [gammaCube,gammaPassRateCell] = matRad_gammaIndex(cube1,cube2,resolution,criteria,slice,n,localglobal,cst)
%   [gammaCube,gammaPassRateCell,gammaFig] = matRad_gammaIndex(cube1,cube2,resolution,criteria,slice,n,localglobal,cst,ct)
%
% input
%   cube1:         dose cube as an M x N x O array
%   cube2:         dose cube as an M x N x O array
%   resolution:    resolution of the cubes [mm/voxel]
%   criteria:      [1x2] vector specifying the distance to agreement
%                  criterion; first element is percentage difference,
%                  second element is distance [mm]
%   slice:         (optional) slice in cube1/2 that will be visualized
%   n:             (optional) number of interpolations. there will be 2^n-1
%                  interpolation points. The maximum suggested value is 3.
%   localglobal:   (optional) parameter to choose between 'global' and 'local'
%                  normalization
%   cst:           list of interessing volumes inside the patient
%   ct:            (optional) matRad ct struct used for slice visualization
%
% output
%
%   gammaCube:          result of gamma index calculation
%   gammaPassRateCell:  rate of voxels passing the specified gamma criterion
%                  as scalar whole-cube rate when 'cst' is omitted, or
%                  evaluated for every structure listed in 'cst'.
%                  note that only voxels exceeding the dose threshold are
%                  considered.
%   gammaFig:           handle to gamma visualization figure if slice was provided
%
% References
%   [1]  http://www.ncbi.nlm.nih.gov/pubmed/9608475
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[env, ~] = matRad_getEnvironment();
gammaFig = [];

if ~exist('slice','var')
    slice = [];
end

if isCstInput(slice) && ~exist('cst','var')
    cst = slice;
    slice = [];
end

if isGammaNormalizationInput(slice) && ~exist('localglobal','var')
    localglobal = char(slice);
    slice = [];
end

if exist('n','var') && isCstInput(n) && ~exist('cst','var')
    cst = n;
    n = [];
end

if exist('localglobal','var') && isCstInput(localglobal) && ~exist('cst','var')
    cst = localglobal;
    localglobal = [];
end

if exist('cst','var') && isstruct(cst) && ~exist('ct','var')
    ct = cst;
    cst = [];
end

if ~exist('cst','var')
    cst = [];
end

if ~exist('ct','var')
    ct = [];
end

% set parameters for gamma index calculation
if exist('criteria','var')
    relDoseThreshold = criteria(1); % in [%]
    dist2AgreeMm     = criteria(2); % in [mm]
else
    dist2AgreeMm     = 3; % in [mm]
    relDoseThreshold = 3; % in [%]
end

% set parameters for gamma index calculation
if ~exist('n','var') || isempty(n)
    n = 0;
end

% set global or local gamma evaluation
if ~exist('localglobal','var') || isempty(localglobal)
    localglobal = 'global';
end

% check if cubes consistent
if ~isequal(size(cube1),size(cube2))
   error('dose cubes must be the same size\n');
end

% define search nighborhood and resolution
if round(n) ~= n
    error('n must be an integer value');
else
    resolution = resolution./(2.^n);
end

neighborhood = 1.5* dist2AgreeMm; % [mm]
searchX = ceil(neighborhood./resolution(1));
searchY = ceil(neighborhood./resolution(2));
searchZ = ceil(neighborhood./resolution(3));

% init cube
% cut all the zeros surrounding the interesting part
cut1 = any(any(cube1,2),3);
cut2 = any(any(cube1,1),3);
cut3 = any(any(cube1,1),2);

% avoids that "zero-slides" are deleted between two interest regions
k = find(cut1); cut1(k(1):k(end)) = 1;
k = find(cut2); cut2(k(1):k(end)) = 1;
k = find(cut2); cut2(k(1):k(end)) = 1;

% copy cubes
cubex1c = cube1;
cubex2c = cube2;

% cut cubes
cubex1c( ~cut1, :, :) = []; % rows
cubex1c( :, ~cut2, :) = []; % columns
cubex1c( :, :, ~cut3) = []; % slices
cubex2c( ~cut1, :, :) = []; % rows
cubex2c( :, ~cut2, :) = []; % columns
cubex2c( :, :, ~cut3) = []; % slices

% pad cubes
cubex2 = zeros([size(cubex2c,1)+2*searchX size(cubex2c,2)+2*searchY size(cubex2c,3)+2*searchZ]);
cubex2((1+searchX):(end-searchX), ...
       (1+searchY):(end-searchY), ...
       (1+searchZ):(end-searchZ)) = cubex2c;

% interpolate if necessary
if n > 0
    switch env
        case 'MATLAB'
            cubex2 = interp3(cubex2,n,'cubic');
        case 'OCTAVE'
            cubex2 = interp3(cubex2,n,'linear');
    end
end

% set up temporary cubes required for calculation
gammaCubeSq = inf*ones(size(cubex1c));

% adjust dose threshold
if strcmp(localglobal,'local')
    doseThreshold = cubex1c .* relDoseThreshold/100;
elseif strcmp(localglobal,'global')
    doseThreshold = relDoseThreshold/100 * max(cube1(:));
end

% search for min
for i = -searchX:searchX
    for j = -searchY:searchY
        for k = -searchZ:searchZ

            delta_sq = ((i*resolution(1))^2 + ...
                        (j*resolution(2))^2 + ...
                        (k*resolution(3))^2) / dist2AgreeMm^2;

            tmpCube = cubex1c - cubex2((1+((2^n)*searchX)+i) : 2^n : (end-((2^n)*searchX)+i), ...
                                       (1+((2^n)*searchY)+j) : 2^n : (end-((2^n)*searchY)+j), ...
                                       (1+((2^n)*searchZ)+k) : 2^n : (end-((2^n)*searchZ)+k));

            tmpCube = tmpCube.^2 ./ doseThreshold.^2 + delta_sq;

            gammaCubeSq = min(gammaCubeSq,tmpCube);


        end
    end

%     display '.';

end

% evaluate gamma cube and set to zero all the voxel that contain an
% infinite value
gammaCubeSqx                 = zeros(size(cube1));
gammaCubeSqx(cut1,cut2,cut3) = gammaCubeSq;
gammaCube                    = sqrt(gammaCubeSqx);

% set values where we did not compute gamma keeping inf in the cube to zero
gammaCube(cube1<=0 & cube2<=0) = 0;

% compute gamma pass rate
doseIx          = cube1 > relDoseThreshold/100*max(cube1(:)) | cube2 > relDoseThreshold/100*max(cube2(:));
numOfPassGamma  = sum(gammaCube(doseIx) < 1);
gammaPassRate   = 100 * numOfPassGamma / sum(doseIx(:));
gammaPassRateCell = gammaPassRate;
refScen = getRefScen(ct);

% compute stats for all segmented sturtures
if ~isempty(cst)
    gammaPassRateCell = cell(1,2);
    gammaPassRateCell{1,1} = 'Whole CT';
    gammaPassRateCell{1,2} = gammaPassRate;

    for i = 1:size(cst, 1)
        volume = cst{i,4}{1,1}; % indices of voxels of the interesting volume
        doseIxVol = false(size(doseIx));
        doseIxVol(volume) = doseIx(volume);
        numOfPassGammaVol  = sum(gammaCube(doseIxVol) < 1);
        gammaPassRateVol   = 100 * numOfPassGammaVol / sum(doseIxVol(:));
        if isnan(gammaPassRateVol) % remove organs not receiving any dose (resulting in NaN)
            continue;
        end
        gammaPassRateCell{end+1, 1} = cst{i,2};
        gammaPassRateCell{end, 2} = gammaPassRateVol;

    end
end

% visualize if applicable
if ~isempty(slice)
    gammaFig = figure;
    set(gammaFig,'Color',[1 1 1]);
    myColormap = matRad_getColormap('gammaIndex');
    gammaWindow = [0 2.01];
    [gammaPlotCube,myColormap] = matRad_prepareLimitedIndexPlot(gammaCube,gammaWindow,myColormap);

    if canUseSliceWrapper(ct)
        gammaFig.Position(3:4) = [400 400];
        matRad_plotSliceWrapper(gca,ct,cst,refScen,gammaPlotCube, ...
            3,slice,[],[],colorcube,myColormap,gammaWindow,[],[], ...
            'Gamma Index',[],'LineWidth',1.5);
    else
        imagesc(gammaPlotCube(:,:,slice),gammaWindow)
        colormap(gca,myColormap);
        colorbar
    end

    subtitle({formatGammaSubtitleText(gammaPassRate,relDoseThreshold,dist2AgreeMm); ...
        ['with ' num2str(2^n-1) ' interpolation points']});
    matRad_addFigureTextbox(gammaFig,formatPercentMetricText('Gamma-PR',gammaPassRate), ...
        'tag','matRadMetricTextbox');
end

end

function tf = canUseSliceWrapper(ct)
tf = isstruct(ct) && isfield(ct,'cubeHU') && iscell(ct.cubeHU) && ...
    isfield(ct,'resolution') && (isfield(ct,'cubeDim') || isfield(ct,'dimensions'));
end

function tf = isCstInput(value)
tf = iscell(value);
end

function tf = isGammaNormalizationInput(value)
tf = (ischar(value) || (isstring(value) && isscalar(value))) && ...
    any(strcmpi(char(value),{'global','local'}));
end

function refScen = getRefScen(ct)
if isstruct(ct) && isfield(ct,'refScen')
    refScen = ct.refScen;
else
    refScen = 1;
end
end

function subtitleText = formatGammaSubtitleText(gammaPassRate,relDoseThreshold,dist2AgreeMm)
subtitleText = [num2str(gammaPassRate,5) '% of points pass gamma criterion (' ...
    num2str(relDoseThreshold) '% / ' ...
    num2str(dist2AgreeMm) 'mm)'];
end

function metricText = formatPercentMetricText(label,value)
if isempty(value) || ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    metricText = [label ' = n/a'];
else
    metricText = sprintf('%s = %.2f%%',label,value);
end
end
