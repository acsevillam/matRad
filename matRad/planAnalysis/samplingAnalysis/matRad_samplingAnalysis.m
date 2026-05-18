function [cstStat, doseStat, meta, figures] = ...
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
%                       - 'evaluationMode': 'perFraction' or 'total' for
%                         figure display
%                       - 'meanDoseWindow': 1x2 mean dose figure window
%                       - 'stdDoseWindow': 1x2 standard deviation figure
%                         window
%                       - 'doseDifferenceWindow': 1x2 signed dose
%                         difference figure window
%                       - 'Percentiles':    vector with desired percentiles
%                       between (0,1)
%
%
% output:
%   cstStat         structure-wise statistics (mean, max, percentiles, ...)
%   doseStat        dose-wise statistics (mean, max, percentiles, ...)
%   meta            contains additional information about sampling analysis
%   figures         formatted sampling figures if a slice was requested
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
figures = matRad_initializeSamplingAnalysisFigures();

meta = matRad_parseSamplingAnalysisInput(ct, pln, varargin);

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

quantityVis = matRad_resolveDoseAnalysisQuantity(resultGUInomScen, pln, '');
doseCube = resultGUInomScen.(quantityVis);

analysisContext = matRad_buildSamplingAnalysisContext(ct, pln, caSampRes, mSampDose, ...
                                                      doseStat, quantityVis, vProb);
meta.analysisContext = analysisContext;
doseStat.analysisContext = analysisContext;
doseStat.nominalAnalysis = matRad_buildSamplingDoseCubeAnalysis( ...
                                                                'nominal', doseCube, ...
                                                                ['resultGUInomScen.' quantityVis], ...
                                                                meta.meanDoseWindow, meta, ...
                                                                analysisContext, pln);
doseStat.meanAnalysis = matRad_buildSamplingDoseCubeAnalysis( ...
                                                             'mean', doseStat.meanCubeW, 'doseStat.meanCubeW', ...
                                                             meta.meanDoseWindow, meta, analysisContext, pln);
doseStat.stdAnalysis = matRad_buildSamplingDoseCubeAnalysis( ...
                                                            'std', doseStat.stdCubeW, 'doseStat.stdCubeW', ...
                                                            meta.stdDoseWindow, meta, analysisContext, pln);

doseStat.gammaAnalysis = matRad_runSamplingGammaAnalysis( ...
                                                         doseCube, doseStat, meta, ct, cstEval, matRad_cfg);

[doseStat.robustnessAnalysis, figures.robustness.index1, figures.robustness.index2] = ...
    matRad_runSamplingRobustnessAnalysis(doseStat, doseCube, cstEval, pln, ...
                                         meta, ct);

figures.mean = matRad_plotSamplingDoseCubeAnalysis(doseStat.meanAnalysis, ...
                                                   doseStat.meanCubeW, ct, cstEval, ...
                                                   meta.slice, 'plane', meta.plane);
figures.std = matRad_plotSamplingDoseCubeAnalysis(doseStat.stdAnalysis, ...
                                                  doseStat.stdCubeW, ct, cstEval, ...
                                                  meta.slice, 'plane', meta.plane);
figures.nominal = matRad_plotSamplingDoseCubeAnalysis(doseStat.nominalAnalysis, ...
                                                      doseCube, ct, cstEval, ...
                                                      meta.slice, 'plane', meta.plane);
[doseStat.expectedDoseDifferenceAnalysis, figures.doseDifference] = ...
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

function meta = matRad_parseSamplingAnalysisInput(ct, pln, args)
matRadCfg = MatRad_Config.instance();
p = inputParser;
p.CaseSensitive = false;
p.addParameter('gammaCriterion', [2 2], @matRad_isPositiveTwoVector);
p.addParameter('robustnessCriteria', [5 5], @matRad_isPositiveTwoVector);
p.addParameter('robustnessTargetMode', 'all', @matRad_isScalarText);
p.addParameter('robustnessTargets', [], @matRad_isStructureSelection);
p.addParameter('expectedDoseDifferenceTolerance', 0, @matRad_isNonNegativeScalar);
p.addParameter('evaluationMode', 'perFraction', @matRad_isScalarText);
p.addParameter('meanDoseWindow', [], @matRad_isDoseWindow);
p.addParameter('stdDoseWindow', [], @matRad_isDoseWindow);
p.addParameter('doseDifferenceWindow', [], @matRad_isDoseWindow);
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

[~, meta.displayEvaluationMode, meta.displayScale] = matRad_convertToEvaluationMode( ...
                                                                                    1, pln, meta.evaluationMode);
meta.doseDifferenceWindowBase = matRad_convertFromEvaluationMode( ...
                                                                 meta.doseDifferenceWindow, pln, ...
                                                                 meta.displayEvaluationMode);
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

function tf = matRad_isDoseWindow(value)
tf = isempty(value) || (isnumeric(value) && isvector(value) && numel(value) == 2 && ...
                        all(isfinite(value(:))) && value(2) > value(1));
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

function figures = matRad_initializeSamplingAnalysisFigures()
figures.mean = [];
figures.std = [];
figures.robustness.index1 = [];
figures.robustness.index2 = [];
figures.nominal = [];
figures.doseDifference = [];
end

function analysis = matRad_buildSamplingDoseCubeAnalysis(analysisType, doseCube, ...
                                                         cubeName, displayWindow, meta, ...
                                                         analysisContext, pln)
analysis.status = matRad_resolveSamplingDoseAnalysisStatus(analysisType, ...
                                                           analysisContext.sampleCoverageFraction);
analysis.cubeName = cubeName;
analysis.quantity = analysisContext.quantity;
analysis.evaluationModeBase = analysisContext.evaluationModeBase;
analysis.displayEvaluationMode = meta.displayEvaluationMode;
analysis.displayScale = meta.displayScale;
analysis.displayDoseWindow = matRad_resolveSamplingDisplayDoseWindow( ...
                                                                     doseCube, displayWindow, ...
                                                                     meta.displayScale, analysisType);
analysis.doseWindow = matRad_convertFromEvaluationMode(analysis.displayDoseWindow, pln, ...
                                                       meta.displayEvaluationMode);
analysis.colorBarLabel = matRad_getSamplingDoseColorBarLabel(analysisType, ...
                                                             analysisContext.quantity);
analysis.title = matRad_getSamplingDoseFigureTitle(analysisType);
end

function status = matRad_resolveSamplingDoseAnalysisStatus(analysisType, sampleCoverageFraction)
if strcmp(analysisType, 'nominal') || sampleCoverageFraction == 1
    status = 'computedFullCube';
else
    status = 'computedPartialMask';
end
end

function displayDoseWindow = matRad_resolveSamplingDisplayDoseWindow(doseCube, ...
                                                                     displayWindow, displayScale, ...
                                                                     analysisType)
if ~isempty(displayWindow)
    displayDoseWindow = displayWindow(:)';
    return
end

displayCube = doseCube .* displayScale;
finiteValues = displayCube(isfinite(displayCube));
if isempty(finiteValues)
    displayDoseWindow = [0 1];
    return
end

maxValue = max(finiteValues(:));
switch analysisType
    case 'std'
        minValue = 0;
    otherwise
        minValue = min(0, min(finiteValues(:)));
end

if maxValue <= minValue
    maxValue = minValue + 1;
end
displayDoseWindow = [minValue maxValue];
end

function label = matRad_getSamplingDoseColorBarLabel(analysisType, quantity)
unitLabel = matRad_getSamplingDoseUnitLabel(quantity);
switch analysisType
    case 'nominal'
        label = ['Nominal dose ' unitLabel];
    case 'std'
        label = ['Dose standard deviation ' unitLabel];
    otherwise
        label = ['Expected dose ' unitLabel];
end
end

function unitLabel = matRad_getSamplingDoseUnitLabel(quantity)
if strncmp(char(quantity), 'RBExD', 5) || strncmp(char(quantity), 'RBExDose', 8)
    unitLabel = '[Gy(RBE)]';
else
    unitLabel = '[Gy]';
end
end

function titleText = matRad_getSamplingDoseFigureTitle(analysisType)
switch analysisType
    case 'nominal'
        titleText = 'Nominal dose';
    case 'std'
        titleText = 'Sampled dose standard deviation';
    otherwise
        titleText = 'Sampled expected dose';
end
end

function gammaAnalysis = matRad_runSamplingGammaAnalysis(doseCube, doseStat, meta, ct, cst, matRadCfg)
analysisContext = meta.analysisContext;
quantityVis = analysisContext.quantity;

gammaAnalysis.cube1Name = ['resultGUInomScen.' quantityVis];
gammaAnalysis.cube1 = doseCube;
gammaAnalysis.cube2 = doseStat.meanCubeW;
gammaAnalysis.cube2Name = 'doseStat.meanCubeW';
gammaAnalysis.scope = matRad_resolveSamplingGammaScope(doseStat.sampleCoverageFraction);
gammaAnalysis.status = 'pending';
gammaAnalysis.reason = '';
gammaAnalysis.doseAgreement = meta.gammaCriterion(1);
gammaAnalysis.distAgreement = meta.gammaCriterion(2);
gammaAnalysis = matRad_attachSamplingAnalysisContext(gammaAnalysis, analysisContext, false);

matRadCfg.dispInfo(['matRad: Performing gamma index analysis with parameters ', ...
                    num2str(meta.gammaCriterion), '[%% mm] \n']);
cstGamma = matRad_prepareReferenceCstForGamma(cst, ct);
[gammaCube1, gammaCube2, gammaScopeMask] = matRad_prepareSamplingGammaCubes( ...
                                                                            doseCube, doseStat.meanCubeW, ...
                                                                            doseStat.sampleMask, gammaAnalysis.scope);
[gammaAnalysis.gammaCube, ~] = matRad_gammaIndex( ...
                                                 gammaCube1, gammaCube2, ...
                                                 [ct.resolution.x ct.resolution.y ct.resolution.z], ...
                                                 meta.gammaCriterion, [], 0, 'global', cstGamma);
gammaAnalysis.gammaCube(~gammaScopeMask) = NaN;
gammaPassRateCell = matRad_calcSamplingGammaPassRates( ...
                                                       gammaAnalysis.gammaCube, doseCube, ...
                                                       doseStat.meanCubeW, meta.gammaCriterion, ...
                                                       gammaScopeMask, cstGamma, gammaAnalysis.scope);

gammaAnalysis.status = matRad_getSamplingGammaComputedStatus(gammaAnalysis.scope);
gammaAnalysis.reason = matRad_getSamplingGammaReason(gammaAnalysis.scope);
gammaAnalysis.gammaPassRateCell = gammaPassRateCell;
gammaAnalysis.gammaPassRate = matRad_getPrimaryGammaPassRate(gammaPassRateCell);
end

function scope = matRad_resolveSamplingGammaScope(sampleCoverageFraction)
if sampleCoverageFraction == 1
    scope = 'wholeCt';
else
    scope = 'sampledVoxels';
end
end

function [gammaCube1, gammaCube2, gammaScopeMask] = ...
    matRad_prepareSamplingGammaCubes(doseCube, meanCube, sampleMask, scope)
gammaCube1 = doseCube;
gammaCube2 = meanCube;

if strcmp(scope, 'wholeCt')
    gammaScopeMask = true(size(doseCube));
else
    gammaScopeMask = sampleMask & isfinite(doseCube) & isfinite(meanCube);
end

gammaCube1(~gammaScopeMask | ~isfinite(gammaCube1)) = 0;
gammaCube2(~gammaScopeMask | ~isfinite(gammaCube2)) = 0;
end

function status = matRad_getSamplingGammaComputedStatus(scope)
if strcmp(scope, 'wholeCt')
    status = 'computedFullCube';
else
    status = 'computedSampledVoxels';
end
end

function reason = matRad_getSamplingGammaReason(scope)
if strcmp(scope, 'wholeCt')
    reason = '';
else
    reason = 'Gamma evaluated on sampled voxels because sampling does not cover whole CT.';
end
end

function gammaPassRateCell = matRad_calcSamplingGammaPassRates(gammaCube, cube1, cube2, ...
                                                               gammaCriterion, scopeMask, ...
                                                               cst, scope)
doseIx = matRad_getSamplingGammaDoseMask(cube1, cube2, gammaCriterion(1), scopeMask);
gammaPassRateCell = cell(1, 2);
gammaPassRateCell{1, 1} = matRad_getSamplingGammaScopeLabel(scope);
gammaPassRateCell{1, 2} = matRad_calcSamplingGammaPassRate(gammaCube, doseIx);

for cstIx = 1:size(cst, 1)
    volume = cst{cstIx, 4}{1, 1};
    doseIxVol = false(size(doseIx));
    doseIxVol(volume) = doseIx(volume);
    gammaPassRateVol = matRad_calcSamplingGammaPassRate(gammaCube, doseIxVol);
    if isnan(gammaPassRateVol)
        continue
    end
    gammaPassRateCell{end + 1, 1} = cst{cstIx, 2};
    gammaPassRateCell{end, 2} = gammaPassRateVol;
end
end

function doseIx = matRad_getSamplingGammaDoseMask(cube1, cube2, relDoseThreshold, scopeMask)
cube1Values = cube1(scopeMask & isfinite(cube1));
cube2Values = cube2(scopeMask & isfinite(cube2));
if isempty(cube1Values)
    maxCube1 = 0;
else
    maxCube1 = max(cube1Values(:));
end
if isempty(cube2Values)
    maxCube2 = 0;
else
    maxCube2 = max(cube2Values(:));
end

doseIx = scopeMask & ...
         ((cube1 > relDoseThreshold / 100 * maxCube1) | ...
          (cube2 > relDoseThreshold / 100 * maxCube2));
end

function label = matRad_getSamplingGammaScopeLabel(scope)
if strcmp(scope, 'wholeCt')
    label = 'Whole CT';
else
    label = 'Sampled voxels';
end
end

function gammaPassRate = matRad_calcSamplingGammaPassRate(gammaCube, doseIx)
if ~any(doseIx(:))
    gammaPassRate = NaN;
else
    gammaPassRate = 100 * sum(gammaCube(doseIx) < 1) / sum(doseIx(:));
end
end

function gammaPassRate = matRad_getPrimaryGammaPassRate(gammaPassRateCell)
if iscell(gammaPassRateCell) && ~isempty(gammaPassRateCell)
    gammaPassRate = gammaPassRateCell{1, 2};
else
    gammaPassRate = gammaPassRateCell;
end
end

function [robustnessAnalysis, robustnessFig1, robustnessFig2] = ...
    matRad_runSamplingRobustnessAnalysis(doseStat, doseCube, cst, pln, meta, ct)
analysisContext = meta.analysisContext;
robustnessAnalysis = matRad_samplingRobustnessAnalysis( ...
                                                       doseStat.meanCubeW, doseStat.stdCubeW, ...
                                                       meta.robustnessCriteria, ct, cst, pln, ...
                                                       [], ...
                                                       'robustnessTargetMode', meta.robustnessTargetMode, ...
                                                       'robustnessTargets', meta.robustnessTargets, ...
                                                       'sampleMask', doseStat.sampleMask);
robustnessAnalysis.sourceCubeName = ['resultGUInomScen.' analysisContext.quantity];
robustnessAnalysis.sourceCube = doseCube;
robustnessAnalysis = matRad_attachSamplingAnalysisContext(robustnessAnalysis, analysisContext, false);
robustnessFig1 = [];
robustnessFig2 = [];
if ~isempty(meta.slice)
    robustnessFig1 = matRad_plotSamplingRobustnessAnalysis(robustnessAnalysis, ct, cst, ...
                                                           meta.slice, 'method', 'index1', ...
                                                           'plane', meta.plane);
    robustnessFig2 = matRad_plotSamplingRobustnessAnalysis(robustnessAnalysis, ct, cst, ...
                                                           meta.slice, 'method', 'index2', ...
                                                           'plane', meta.plane);
end
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
                                                                        analysisContext.quantity], ...
                                                                       'doseWindow', ...
                                                                       meta.doseDifferenceWindowBase);
expectedDoseDifferenceAnalysis = matRad_attachSamplingAnalysisContext(expectedDoseDifferenceAnalysis, ...
                                                                      analysisContext, true);
if ~isempty(meta.slice)
    expectedDoseDifferenceFig = ...
        matRad_plotExpectedDoseDifferenceAnalysis(expectedDoseDifferenceAnalysis, ct, cst, ...
                                                  meta.slice, 'plane', meta.plane, ...
                                                  'doseWindow', meta.doseDifferenceWindow, ...
                                                  'displayScale', meta.displayScale);
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
