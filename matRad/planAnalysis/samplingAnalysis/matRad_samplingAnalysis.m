function [cstStat, doseStat, meta, gammaFig, robustnessFig1, robustnessFig2, expectedDoseDifferenceFig] = ...
    matRad_samplingAnalysis(ct, cst, pln, caSampRes, mSampDose, resultGUInomScen, varargin)
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
%   gammaFig        gamma analysis figure if a slice was requested
%   robustnessFig1  robustness index1 figure if a slice was requested
%   robustnessFig2  robustness index2 figure if a slice was requested
%   expectedDoseDifferenceFig  expected dose difference figure if a slice was requested
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
gammaFig = [];
robustnessFig1 = [];
robustnessFig2 = [];
expectedDoseDifferenceFig = [];

meta = matRad_parseSamplingAnalysisInput(ct, varargin);

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

analysisContext = matRad_buildSamplingAnalysisContext(ct, pln, caSampRes, mSampDose, ...
                                                      doseStat, quantityVis, vProb);
meta.analysisContext = analysisContext;
doseStat.analysisContext = analysisContext;

[doseStat.gammaAnalysis, gammaFig] = matRad_runSamplingGammaAnalysis( ...
                                                                     doseCube, doseStat, meta, ct, cstEval, matRad_cfg);

[doseStat.robustnessAnalysis, robustnessFig1, robustnessFig2] = ...
    matRad_runSamplingRobustnessAnalysis(doseStat, doseCube, cstEval, pln, ...
                                         meta, ct);

[doseStat.expectedDoseDifferenceAnalysis, expectedDoseDifferenceFig] = ...
    matRad_runSamplingExpectedDoseDifferenceAnalysis(mSampDose, doseStat, doseCube, ...
                                                     meta, ct, cstEval);

end

function statCheck = matRad_checkSampIntegrity(multScen)
% check integrity of scenario analysis (i.e. check number of scenarios)
if multScen.numScenarios() > 20
    statCheck = true;
else
    statCheck = false;
end
end

function meta = matRad_parseSamplingAnalysisInput(ct, args)
matRadCfg = MatRad_Config.instance();
p = inputParser;
p.CaseSensitive = false;
p.addParameter('gammaCriterion', [2 2], @matRad_isPositiveTwoVector);
p.addParameter('robustnessCriteria', [5 5], @matRad_isPositiveTwoVector);
p.addParameter('robustnessTargetMode', 'all', @matRad_isScalarText);
p.addParameter('robustnessTargets', [], @matRad_isStructureSelection);
p.addParameter('expectedDoseDifferenceTolerance', 0, @matRad_isNonNegativeScalar);
p.addParameter('ctScenProb', [], @matRad_isCtScenarioProbabilityOverride);
p.addParameter('slice', [], @matRad_isOptionalPositiveIntegerScalar);
p.addParameter('plane', 3, @matRad_isValidPlane);
p.addParameter('percentiles', [0.01 0.05 0.125 0.25 0.5 0.75 0.875 0.95 0.99], ...
               @matRad_isPercentileVector);

parse(p, args{:});
meta = p.Results;

if ~isempty(meta.slice) && meta.slice > ct.cubeDim(3)
    matRadCfg.dispError('slice must be between 1 and %d.', ct.cubeDim(3));
end
end

function tf = matRad_isPositiveTwoVector(value)
tf = numel(value) == 2 && isnumeric(value) && all(value > 0);
end

function tf = matRad_isScalarText(value)
tf = ischar(value) || (isstring(value) && isscalar(value));
end

function tf = matRad_isStructureSelection(value)
tf = isempty(value) || isnumeric(value) || ischar(value) || isstring(value) || iscell(value);
end

function tf = matRad_isNonNegativeScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0;
end

function tf = matRad_isCtScenarioProbabilityOverride(value)
tf = isempty(value) || (isnumeric(value) && ismatrix(value) && size(value, 2) == 2 && ...
                        all(isfinite(value(:))) && all(value(:, 2) >= 0));
end

function tf = matRad_isOptionalPositiveIntegerScalar(value)
tf = isempty(value) || (isnumeric(value) && isscalar(value) && isfinite(value) && ...
                        round(value) == value && value >= 1);
end

function tf = matRad_isValidPlane(value)
tf = isnumeric(value) && isscalar(value) && any(value == [1 2 3]);
end

function tf = matRad_isPercentileVector(value)
tf = (isscalar(value) || isvector(value)) && isnumeric(value) && all(value > 0 & value < 1);
end

function [gammaAnalysis, gammaFig] = matRad_runSamplingGammaAnalysis(doseCube, doseStat, meta, ct, cst, matRadCfg)
gammaFig = [];
analysisContext = meta.analysisContext;
quantityVis = analysisContext.quantity;

gammaAnalysis.cube1Name = ['resultGUInomScen.' quantityVis];
gammaAnalysis.cube1 = doseCube;
gammaAnalysis.cube2 = doseStat.meanCubeW;
gammaAnalysis.cube2Name = 'doseStat.meanCubeW';
gammaAnalysis.scope = 'wholeCt';
gammaAnalysis.status = 'pending';
gammaAnalysis.reason = '';
gammaAnalysis.doseAgreement = meta.gammaCriterion(1);
gammaAnalysis.distAgreement = meta.gammaCriterion(2);
gammaAnalysis = matRad_attachSamplingAnalysisContext(gammaAnalysis, analysisContext, false);

if doseStat.sampleCoverageFraction < 1
    gammaAnalysis = matRad_markGammaSkippedForPartialSampling(gammaAnalysis, ct, doseStat, ...
                                                              matRadCfg);
    return
end

matRadCfg.dispInfo(['matRad: Performing gamma index analysis with parameters ', ...
                    num2str(meta.gammaCriterion), '[%% mm] \n']);
cstGamma = matRad_prepareReferenceCstForGamma(cst, ct);
[gammaAnalysis.gammaCube, gammaPassRateCell] = matRad_gammaIndex( ...
                                                                 doseCube, doseStat.meanCubeW, ...
                                                                 [ct.resolution.x ct.resolution.y ct.resolution.z], ...
                                                                 meta.gammaCriterion, meta.slice, 0, 'global', cstGamma);

gammaAnalysis.status = 'computedFullCube';
gammaAnalysis.gammaPassRateCell = gammaPassRateCell;
gammaAnalysis.gammaPassRate = matRad_getWholeCtGammaPassRate(gammaPassRateCell);

if ~isempty(meta.slice)
    gammaFig = gcf;
end
end

function gammaAnalysis = matRad_markGammaSkippedForPartialSampling(gammaAnalysis, ct, doseStat, matRadCfg)
gammaAnalysis.status = 'skippedPartialSampling';
gammaAnalysis.reason = 'Gamma whole-CT analysis requires sampled statistics for every CT voxel.';
gammaAnalysis.gammaCube = NaN(ct.cubeDim);
gammaAnalysis.gammaPassRate = [];
gammaAnalysis.gammaPassRateCell = {};
matRadCfg.dispWarning(['Skipping gamma index analysis because sampling only ', ...
                       'covers %.2f%% of CT voxels.\n'], ...
                      100 * doseStat.sampleCoverageFraction);
end

function gammaPassRate = matRad_getWholeCtGammaPassRate(gammaPassRateCell)
if iscell(gammaPassRateCell) && ~isempty(gammaPassRateCell)
    gammaPassRate = gammaPassRateCell{1, 2};
else
    gammaPassRate = gammaPassRateCell;
end
end

function [robustnessAnalysis, robustnessFig1, robustnessFig2] = ...
    matRad_runSamplingRobustnessAnalysis(doseStat, doseCube, cst, pln, meta, ct)
analysisContext = meta.analysisContext;
[robustnessAnalysis, robustnessFig1, robustnessFig2] = ...
    matRad_samplingRobustnessAnalysis(doseStat.meanCubeW, doseStat.stdCubeW, ...
                                      meta.robustnessCriteria, ct, cst, pln, ...
                                      meta.slice, ...
                                      'robustnessTargetMode', meta.robustnessTargetMode, ...
                                      'robustnessTargets', meta.robustnessTargets, ...
                                      'sampleMask', doseStat.sampleMask);
robustnessAnalysis.sourceCubeName = ['resultGUInomScen.' analysisContext.quantity];
robustnessAnalysis.sourceCube = doseCube;
robustnessAnalysis = matRad_attachSamplingAnalysisContext(robustnessAnalysis, analysisContext, false);
end

function [expectedDoseDifferenceAnalysis, expectedDoseDifferenceFig] = ...
    matRad_runSamplingExpectedDoseDifferenceAnalysis(mSampDose, doseStat, doseCube, meta, ct, cst)
expectedDoseDifferenceFig = [];
analysisContext = meta.analysisContext;
expectedDoseDifferenceAnalysis = matRad_expectedDoseDifferenceAnalysis( ...
                                                                       mSampDose, doseCube, ...
                                                                       'sampledVoxelIndices', ...
                                                                       analysisContext.sampledVoxelIndices, ...
                                                                       'sampleMask', analysisContext.sampleMask, ...
                                                                       'scenWeights', analysisContext.scenWeights, ...
                                                                       'tolerance', ...
                                                                       meta.expectedDoseDifferenceTolerance, ...
                                                                       'quantity', analysisContext.quantity, ...
                                                                       'evaluationModeBase', ...
                                                                       analysisContext.evaluationModeBase, ...
                                                                       'referenceName', ...
                                                                       ['resultGUInomScen.' ...
                                                                        analysisContext.quantity]);
expectedDoseDifferenceAnalysis = matRad_attachSamplingAnalysisContext(expectedDoseDifferenceAnalysis, ...
                                                                      analysisContext, true);
if ~isempty(meta.slice)
    expectedDoseDifferenceFig = ...
        matRad_plotExpectedDoseDifferenceAnalysis(expectedDoseDifferenceAnalysis, ct, cst, ...
                                                  meta.slice, 'plane', meta.plane);
end
end

function analysis = matRad_attachSamplingAnalysisContext(analysis, analysisContext, includeScenarioData)
analysis.analysisGrid = analysisContext.analysisGrid;
analysis.referenceCtScen = analysisContext.referenceCtScen;
analysis.quantity = analysisContext.quantity;
analysis.evaluationModeBase = analysisContext.evaluationModeBase;
analysis.sampleCoverageFraction = analysisContext.sampleCoverageFraction;
analysis.numSamples = analysisContext.numSamples;

if includeScenarioData
    analysis.scenarioIds = analysisContext.scenarioIds;
    analysis.ctScenIds = analysisContext.ctScenIds;
    analysis.doseMapping = analysisContext.doseMapping;
end
end

function cstGamma = matRad_prepareReferenceCstForGamma(cst, ct)
cstGamma = cst;
if isfield(ct, 'refScen')
    refScen = ct.refScen;
else
    refScen = 1;
end

for cstIx = 1:size(cstGamma, 1)
    if numel(cstGamma{cstIx, 4}) >= refScen
        cstGamma{cstIx, 4}{1} = cstGamma{cstIx, 4}{refScen};
    end
end
end
