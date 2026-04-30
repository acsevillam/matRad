function matRad_createAnimationForLatexReport(confidenceValue, ct, cst, slice, meanCube, mRealizations, scenWeights, subIx, outpath, legendColorbar,varargin)
% matRad function to create figures for a GIF animation
%
% call
%   matRad_createAnimationForLatexReport(confidenceValue, ct, cst, slice, ...
%           meanCube, mRealizations, scenWeights, subIx, outpath, legendColorbar)
%
% input
%   confidenceValue: confidence used for visualization
%   ct:              matRad ct struct
%   cst:             matRad cst struct
%   slice:           slice of the ct used for visualization
%   meanCube:        cube holding the mean dose
%   mRealizations:   samples
%   scenWeights:     linear vector of all weights for the individual
%                    scenarios
%   subIx:           voxel indices that are considered during analysis
%   outpath:         output path for files
%   legendColorbar:  colorbar used for the legend
%
% input (optional Name-Value pairs)
%   varargin:        optional Name-Value pairs
%   PrescribedDose:  prescription (per fraction)
%   FramesPerSecond: frames per second for the animation (default 24)
%   Period:          total period [s] for the animation (default 5)
%   FilePrefix:      file prefix (default 'anim')
%
% output
%   -
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

p.addParameter('PrescribedDose',[],@(x) isscalar(x) && isnumeric(x));
p.addParameter('FramesPerSecond',24,@(x) isscalar(x) && isnumeric(x));
p.addParameter('Period',5,@(x) isscalar(x) && isnumeric(x));
p.addParameter('FilePrefix','anim',@ischar);

p.parse(varargin{:});

fps = p.Results.FramesPerSecond;
period = p.Results.Period;
dPres = p.Results.PrescribedDose;

ctDim = size(meanCube);
doseSlice = meanCube(:,:,slice);
doseSliceIx = find(isfinite(doseSlice) & doseSlice ~= 0) + (slice-1)*prod(ctDim(1:2));

if isempty(doseSliceIx)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispWarning('Skipping sampling animation because the selected slice has no finite sampled dose values.\n');
    return;
end

% check whether expDoseSliceIx is subset of subIx
assert(all(ismember(doseSliceIx, subIx)))
assert(issorted(doseSliceIx) && issorted(subIx))

if isempty(dPres)
    finiteDoseSlice = doseSlice(isfinite(doseSlice));
    dPres = 0.95 * max(finiteDoseSlice(:));
end

[~,idxIntoSubIx] = intersect(subIx, doseSliceIx);
mRealizationsSub = mRealizations(idxIntoSubIx,:);

selectIx = doseSliceIx;

scenWeights = scenWeights(:);
scenWeights(~isfinite(scenWeights) | scenWeights < 0) = 0;
if sum(scenWeights) <= 0
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Scenario weights must contain at least one positive finite value.');
end
scenWeights = scenWeights./sum(scenWeights);

[~,nSamples] = size(mRealizationsSub);
if numel(scenWeights) ~= nSamples
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Number of scenario weights does not match number of sampled scenarios.');
end

weightedMean = scenWeights' * mRealizationsSub';
centeredRealizations = bsxfun(@minus,mRealizationsSub',weightedMean);
wCovMat = centeredRealizations' * bsxfun(@times,centeredRealizations,scenWeights);
wCovMat = 0.5 * (wCovMat + wCovMat');

mu = meanCube(selectIx);

%Compute starting vector from confidence value
xr = sqrt(gammaincinv(confidenceValue,numel(mu)/2)*2);


fname = 'anim';
nFrames = period*fps;
samples = matRad_getGaussianOrbitSamples(mu,wCovMat,nFrames,xr);
samplesMax = max(samples(:));
hfAnim = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
%outfile = [fname '.gif'];
alpha = 0.5;
close
figure
set(gcf,'color','w');
for f=1:nFrames
    sampleCube = NaN(size(meanCube));
    sampleCube(selectIx) = samples(:,f);
    doseWindow = [0.01*dPres dPres*1.3];
    isoLevels = [0.1 0.25 0.6 0.9 0.95 1 1.05 1.25]'*dPres;
    matRad_plotSliceWrapper(gca,ct,cst,1,sampleCube,3,slice,0,alpha, ...
        colorcube,jet,doseWindow,isoLevels,[],legendColorbar,false);
    F(f) = getframe(gcf);
    im = frame2im(F(f));
    [imind,cm] = rgb2ind(im,256);

    outfile = fullfile(outpath,[fname '_' num2str(f) '.png']);
    imwrite(imind,cm,outfile,'png');
    delete(gca);
end
close
end
