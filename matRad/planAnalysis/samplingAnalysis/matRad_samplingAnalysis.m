function [cstStat, doseStat, meta] = matRad_samplingAnalysis(ct, cst, pln, caSampRes, mSampDose, resultGUInomScen, varargin)
% matRad uncertainty sampling analysis function
%
% call:
%   [structureStat, doseStat] = samplingAnalysis(ct,cst,subIx,mSampDose,w)
%
% input:
%   ct:                 ct cube
%   cst:                matRad cst struct
%   pln:                matRad's pln struct
%   caSampRes:          cell array of sampling results depicting plan
%                       parameter
%   mSampDose:          matrix holding the sampled doses, each row
%                       corresponds to one dose sample
%   resultGUInomScen:   resultGUI struct of the nominal plan
%   varargin:           optional Name/Value pairs for additional custom
%                       settings
%                       - 'GammaCriterion': 1x2 vector [%  mm]
%                       - 'Percentiles':    vector with desired percentiles
%                       between (0,1)
%
%
% output:
%   cstStat         structure-wise statistics (mean, max, percentiles, ...)
%   doseStat        dose-wise statistics (mean, max, percentiles, ...)
%   meta            contains additional information about sampling analysis
%
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

p = inputParser;
p.addParameter('gammaCriterion', [2 2], @(g) numel(g) == 2 && isnumeric(g) && all(g > 0));
p.addParameter('ctScenProb', [], @(ctScenProb) isempty(ctScenProb) || ...
               (isnumeric(ctScenProb) && ismatrix(ctScenProb) && size(ctScenProb, 2) == 2 && ...
                all(isfinite(ctScenProb(:))) && all(ctScenProb(:, 2) >= 0)));
p.addParameter('percentiles', [0.01 0.05 0.125 0.25 0.5 0.75 0.875 0.95 0.99], ...
               @(x) (isscalar(x) || isvector(x)) && isnumeric(x) && all(x > 0 & x < 1));

parse(p, varargin{:});

meta = p.Results;

meta.sufficientStatistics = matRad_checkSampIntegrity(pln.multScen);
meta.scenWeights = matRad_getSamplingScenarioWeights(pln, numel(caSampRes), meta.ctScenProb);
meta.evaluationModeBase = 'perFraction';

vProb = meta.scenWeights;

%% generate structurewise statistics
if isfield(resultGUInomScen, 'cst') && ~isempty(resultGUInomScen.cst)
    cstEval = resultGUInomScen.cst;
else
    cstEval = cst;
end
cstStat = matRad_calcSamplingStructureStatistics(cstEval, caSampRes, vProb, meta.percentiles);

%% calculate mean and std cube
doseStat = matRad_calcSamplingDoseStatistics(ct, pln, mSampDose, vProb);

% gamma cube
quantityVis = matRad_resolveDoseAnalysisQuantity(resultGUInomScen, pln, '');

doseCube = resultGUInomScen.(quantityVis);

doseStat.gammaAnalysis.cube1Name = ['resultGUInomScen.' quantityVis];

doseStat.gammaAnalysis.cube1 = doseCube;
doseStat.gammaAnalysis.cube2 = doseStat.meanCubeW;
doseStat.gammaAnalysis.cube2Name = 'doseStat.meanCubeW';

matRad_cfg.dispInfo(['matRad: Performing gamma index analysis with parameters', num2str(meta.gammaCriterion), ...
                     '[%% mm] \n']);
doseStat.gammaAnalysis.doseAgreement = meta.gammaCriterion(1);
doseStat.gammaAnalysis.distAgreement = meta.gammaCriterion(2);

doseCubeForGamma = doseStat.meanCubeW;
doseCubeForGamma(~doseStat.sampleMask) = 0;
doseStat.gammaAnalysis.gammaCube = matRad_gammaIndex(doseCube, doseCubeForGamma, ...
                                                     [ct.resolution.x ct.resolution.y ct.resolution.z], ...
                                                     meta.gammaCriterion);

end

function statCheck = matRad_checkSampIntegrity(multScen)
% check integrity of scenario analysis (i.e. check number of scenarios)
if multScen.numScenarios() > 20
    statCheck = true;
else
    statCheck = false;
end
end
