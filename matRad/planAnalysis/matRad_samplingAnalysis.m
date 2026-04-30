function [cstStat, doseStat, meta, gammaFig, robustnessFig1, robustnessFig2] = ...
    matRad_samplingAnalysis(ct,cst,pln,caSampRes,mSampDose,resultGUInomScen,varargin)
% matRad uncertainty sampling analysis function
%
% call
%   [cstStat,doseStat] = matRad_samplingAnalysis(ct,cst,pln,caSampRes,mSampDose,resultGUInomScen)
%   [cstStat,doseStat,meta,gammaFig,robustnessFig1,robustnessFig2] = matRad_samplingAnalysis(ct,cst,pln,caSampRes,mSampDose,resultGUInomScen,...)
%   [...] = matRad_samplingAnalysis(...,'ctScenProb',ctScenProb,'slice',slice)
%
% input
%   ct:                 ct cube
%   cst:                matRad cst struct
%   pln:                matRad's pln struct
%   caSampRes:          cell array of sampling results depicting plan
%                       parameter
%   mSampDose:          sampled dose matrix with one row per voxel and one
%                       column per sampled scenario
%   resultGUInomScen:   resultGUI struct of the nominal plan
%   varargin:           optional Name/Value pairs for additional custom
%                       settings
%                       - 'ctScenProb':     CT scenario probability override
%                                           [ctScenIndex probability]
%                       - 'gammaCriterion': 1x2 vector [%  mm]
%                       - 'robustnessCriteria': 1x2 vector [%  %]
%                       - 'robustnessTargetMode': 'all', 'include', or
%                                           'exclude'
%                       - 'robustnessTargets': target names or cst row
%                                           indices for include/exclude mode
%                       - 'slice':          CT slice used to create figures
%                       Dose statistics are evaluated per fraction.
%                       - 'percentiles':    vector with desired percentiles
%                                           between (0,1)
%
%
% output
%   cstStat:           structure-wise statistics (mean, max, percentiles, ...)
%   doseStat:          dose-wise statistics (mean, max, percentiles, ...)
%   meta:              contains additional information about sampling analysis
%   gammaFig:          gamma analysis figure if slice was provided
%   robustnessFig1:    robustness analysis figure for combined deviation/std
%                      metric if slice was provided
%   robustnessFig2:    robustness analysis figure for binary robustness
%                      criterion if slice was provided
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

%% check integrity of statistics
matRad_cfg = MatRad_Config.instance();
robustnessFig1 = [];
robustnessFig2 = [];

p = inputParser;
p.CaseSensitive = false;
p.addParameter('gammaCriterion',[2 2],@(g) numel(g) == 2 && isnumeric(g) && all(g > 0));
p.addParameter('robustnessCriteria',[5 5],@(r) numel(r) == 2 && isnumeric(r) && all(r > 0));
p.addParameter('robustnessTargetMode','all',@(m) ischar(m) || (isstring(m) && isscalar(m)));
p.addParameter('robustnessTargets',[],@(s) isempty(s) || isnumeric(s) || ischar(s) || isstring(s) || iscell(s));
p.addParameter('ctScenProb',[],@(ctScenProb) isempty(ctScenProb) || ...
    (isnumeric(ctScenProb) && ismatrix(ctScenProb) && size(ctScenProb,2) == 2 && ...
    all(isfinite(ctScenProb(:))) && all(ctScenProb(:,2) >= 0)));
p.addParameter('slice',[],@(slice) isempty(slice) || ...
    (isnumeric(slice) && isscalar(slice) && isfinite(slice) && round(slice) == slice && slice >= 1));
p.addParameter('percentiles',[0.01 0.05 0.125 0.25 0.5 0.75 0.875 0.95 0.99],@(p) (isscalar(p) || isvector(p)) && isnumeric(p) && all(p > 0 & p < 1));

parse(p,varargin{:});

meta = p.Results;
if ~isempty(meta.slice) && meta.slice > ct.cubeDim(3)
    matRad_cfg.dispError('slice must be between 1 and %d.',ct.cubeDim(3));
end

meta.evaluationModeBase = 'perFraction';

meta.sufficientStatistics = matRad_checkSampIntegrity(pln.multScen);
calcRobustness = nargout > 4 || ~any(strcmp(p.UsingDefaults,'robustnessCriteria'));

scenWeights = matRad_getSamplingScenarioWeights(pln,numel(caSampRes),meta.ctScenProb);
meta.scenWeights = scenWeights;

%% generate structurewise statistics
cstStat = matRad_calcSamplingStructureStatistics(cst,caSampRes,scenWeights,meta.percentiles);

%% calculate mean and std cube
doseStat = matRad_calcSamplingDoseStatistics(ct,pln,mSampDose,scenWeights);

% gamma cube
doseCube = resultGUInomScen.(pln.bioParam.quantityVis);
if strncmp(pln.bioParam.quantityVis,'RBExD', 5)
    doseStat.gammaAnalysis.cube1Name = 'resultGUInomScen.RBExD';
else
    doseStat.gammaAnalysis.cube1Name = 'resultGUInomScen.physicalDose';
end

doseStat.gammaAnalysis.cube1 = doseCube;
doseStat.gammaAnalysis.cube2 = doseStat.meanCubeW;
doseStat.gammaAnalysis.cube2Name = 'doseStat.meanCubeW';
doseStat.gammaAnalysis.scope = 'wholeCt';
doseStat.gammaAnalysis.status = 'pending';
doseStat.gammaAnalysis.reason = '';

doseStat.gammaAnalysis.doseAgreement = meta.gammaCriterion(1);
doseStat.gammaAnalysis.distAgreement = meta.gammaCriterion(2);

if doseStat.sampleCoverageFraction < 1
    doseStat.gammaAnalysis.status = 'skippedPartialSampling';
    doseStat.gammaAnalysis.reason = ...
        'Gamma whole-CT analysis requires sampled statistics for every CT voxel.';
    doseStat.gammaAnalysis.gammaCube = NaN(ct.cubeDim);
    doseStat.gammaAnalysis.gammaPassRate = [];
    doseStat.gammaAnalysis.gammaPassRateCell = {};
    gammaFig = [];
    matRad_cfg.dispWarning(['Skipping gamma index analysis because sampling only ', ...
        'covers %.2f%% of CT voxels.\n'],100*doseStat.sampleCoverageFraction);
else
    doseStat.gammaAnalysis.status = 'computedFullCube';
    matRad_cfg.dispInfo(['matRad: Performing gamma index analysis with parameters ', ...
        num2str(meta.gammaCriterion), ' [%% mm] \n']);
    [doseStat.gammaAnalysis.gammaCube,gammaPassRate,gammaFig] = ...
        matRad_gammaIndex(doseCube,doseStat.meanCubeW, ...
        [ct.resolution.x ct.resolution.y ct.resolution.z],meta.gammaCriterion,meta.slice,[],[],cst,ct);

    if iscell(gammaPassRate)
        doseStat.gammaAnalysis.gammaPassRate = gammaPassRate{1,2};
        doseStat.gammaAnalysis.gammaPassRateCell = gammaPassRate;
    else
        doseStat.gammaAnalysis.gammaPassRate = gammaPassRate;
    end
end

if calcRobustness
    if nargout > 5
        [doseStat.robustnessAnalysis,robustnessFig1,robustnessFig2] = ...
            matRad_samplingRobustnessAnalysis(doseStat.meanCubeW,doseStat.stdCubeW, ...
            meta.robustnessCriteria,ct,cst,pln,meta.slice, ...
            'robustnessTargetMode',meta.robustnessTargetMode, ...
            'robustnessTargets',meta.robustnessTargets, ...
            'sampleMask',doseStat.sampleMask);
    elseif nargout > 4
        [doseStat.robustnessAnalysis,robustnessFig1] = ...
            matRad_samplingRobustnessAnalysis(doseStat.meanCubeW,doseStat.stdCubeW, ...
            meta.robustnessCriteria,ct,cst,pln,meta.slice, ...
            'robustnessTargetMode',meta.robustnessTargetMode, ...
            'robustnessTargets',meta.robustnessTargets, ...
            'sampleMask',doseStat.sampleMask);
    else
        doseStat.robustnessAnalysis = ...
            matRad_samplingRobustnessAnalysis(doseStat.meanCubeW,doseStat.stdCubeW, ...
            meta.robustnessCriteria,ct,cst,pln,meta.slice, ...
            'robustnessTargetMode',meta.robustnessTargetMode, ...
            'robustnessTargets',meta.robustnessTargets, ...
            'sampleMask',doseStat.sampleMask);
    end

    if strncmp(pln.bioParam.quantityVis,'RBExD', 5)
        doseStat.robustnessAnalysis.sourceCubeName = 'resultGUInomScen.RBExD';
    else
        doseStat.robustnessAnalysis.sourceCubeName = 'resultGUInomScen.physicalDose';
    end
    doseStat.robustnessAnalysis.sourceCube = doseCube;
end

    % check integrity of scenario analysis (i.e. check number of scenarios)
    function statCheck = matRad_checkSampIntegrity(multScen)

        if multScen.totNumScen > 20
            totalNum = true;
        else
            totalNum = false;
        end

        statCheck = totalNum; % * .... *

    end

end
