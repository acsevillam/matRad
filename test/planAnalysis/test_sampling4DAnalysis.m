function test_suite = test_sampling4DAnalysis

test_functions = localfunctions();

initTestSuite;

function test_samplingDoseStatisticsUseWeightedPopulationStd
ct.cubeDim = [2 2 1];
pln.subIx = [1; 4];
mSampDose = [1 3 5; ...
             2 4 8];
weights = [0.2; 0.3; 0.5];

doseStat = matRad_calcSamplingDoseStatistics(ct, pln, mSampDose, weights);

expectedMean = mean(mSampDose, 2);
expectedStd = sqrt(mean((mSampDose - expectedMean).^2, 2));
expectedMeanW = mSampDose * weights;
expectedStdW = sqrt(((mSampDose - expectedMeanW).^2) * weights);

assertElementsAlmostEqual(doseStat.meanCube(pln.subIx), expectedMean, 'absolute', 1e-12);
assertElementsAlmostEqual(doseStat.stdCube(pln.subIx), expectedStd, 'absolute', 1e-12);
assertElementsAlmostEqual(doseStat.meanCubeW(pln.subIx), expectedMeanW, 'absolute', 1e-12);
assertElementsAlmostEqual(doseStat.stdCubeW(pln.subIx), expectedStdW, 'absolute', 1e-12);
assertEqual(doseStat.sampledVoxelIndices, pln.subIx);
assertElementsAlmostEqual(doseStat.sampleCoverageFraction, 0.5, 'absolute', 1e-12);
assertTrue(doseStat.sampleMask(1));
assertTrue(doseStat.sampleMask(4));
assertTrue(isnan(doseStat.meanCube(2)));
assertTrue(isnan(doseStat.stdCubeW(3)));

function test_samplingStructureStatisticsUseWeightedDvhStd
cst = helper_createCst();
caSampRes = helper_createSamplingResults();
weights = [0.5; 0.25; 0.25];
percentiles = [0.5];

cstStat = matRad_calcSamplingStructureStatistics(cst, caSampRes, weights, percentiles);

dvhMat = zeros(numel(caSampRes), numel(caSampRes(1).dvh(1).volumePoints));
for i = 1:numel(caSampRes)
    dvhMat(i, :) = caSampRes(i).dvh(1).volumePoints;
end
expectedMean = weights(:)' * dvhMat;
expectedStd = sqrt(weights(:)' * (dvhMat - expectedMean).^2);

assertElementsAlmostEqual(cstStat(1).dvhStat.mean.volumePoints, expectedMean, 'absolute', 1e-12);
assertElementsAlmostEqual(cstStat(1).dvhStat.std.volumePoints, expectedStd, 'absolute', 1e-12);
assertElementsAlmostEqual(cstStat(1).w(:), weights, 'absolute', 1e-12);

function test_samplingMappingMetadataDescribesFourDReference
ct = helper_createCtWithPullDvf();
model = matRad_NominalScenario(ct);
doseMapping = matRad_resolveSamplingDoseMapping(ct, model, ct.refScen);

assertTrue(doseMapping.enabled);
assertEqual(doseMapping.refScen, 1);
assertEqual(doseMapping.method, 'pullDirectDoseMapping');

sample = helper_createSamplingResult(2, 2, 1, true);

assertEqual(sample.ctScenId, 2);
assertEqual(sample.refScen, 1);
assertTrue(sample.doseMapping.mapped);
assertEqual(sample.doseMapping.sourceCtScenId, 2);
assertEqual(sample.doseMapping.referenceCtScenId, 1);

function test_samplingAnalysisContextDescribesReferenceGrid
ct = helper_createFourDAnalysisCt();
cst = helper_createCst();
pln = helper_createAnalysisPlan(ct, [1; 2; 1]);
caSampRes = helper_createFourDSamplingResults();
mSampDose = [1 2 4; ...
             2 3 5];
resultGUInomScen = helper_createNominalResult(ct, cst);

[~, doseStat, meta] = matRad_samplingAnalysis(ct, cst, pln, caSampRes, ...
                                              mSampDose, resultGUInomScen);

context = meta.analysisContext;
assertEqual(context.analysisGrid, 'referenceCt');
assertEqual(context.referenceCtScen, 1);
assertEqual(context.quantity, 'physicalDose');
assertEqual(context.evaluationModeBase, 'perFraction');
assertEqual(context.scenarioIds(:), [1; 2; 3]);
assertEqual(context.ctScenIds(:), [1; 2; 1]);
assertElementsAlmostEqual(context.scenWeights(:), [0.2; 0.3; 0.5], 'absolute', 1e-12);
assertEqual(context.sampledVoxelIndices(:), pln.subIx(:));
assertElementsAlmostEqual(context.sampleCoverageFraction, 0.25, 'absolute', 1e-12);

assertEqual(doseStat.analysisContext.analysisGrid, 'referenceCt');
assertEqual(doseStat.gammaAnalysis.analysisGrid, 'referenceCt');
assertEqual(doseStat.robustnessAnalysis.referenceCtScen, 1);
assertEqual(doseStat.expectedDoseDifferenceAnalysis.analysisGrid, 'referenceCt');
assertEqual(doseStat.expectedDoseDifferenceAnalysis.ctScenIds(:), [1; 2; 1]);
assertTrue(doseStat.expectedDoseDifferenceAnalysis.doseMapping(2).mapped);

function test_samplingAnalysisRejectsUnmappedFourDScenario
ct = helper_createFourDAnalysisCt();
cst = helper_createCst();
pln = helper_createAnalysisPlan(ct, [1; 2; 1]);
caSampRes = helper_createFourDSamplingResults();
caSampRes(2).doseMapping.mapped = false;
mSampDose = [1 2 4; ...
             2 3 5];
resultGUInomScen = helper_createNominalResult(ct, cst);

assertExceptionThrown(@() matRad_samplingAnalysis(ct, cst, pln, caSampRes, ...
                                                  mSampDose, resultGUInomScen), ...
                      'matRad:Error');

function test_samplingAnalysisRejectsQuantityMismatch
ct = helper_createAnalysisCt();
cst = helper_createCst();
pln = helper_createAnalysisPlan(ct);
caSampRes = helper_createSamplingResults();
caSampRes(1).analysisQuantity = 'RBExD';
mSampDose = [1 2 4; ...
             2 3 5];
resultGUInomScen = helper_createNominalResult(ct, cst);

assertExceptionThrown(@() matRad_samplingAnalysis(ct, cst, pln, caSampRes, ...
                                                  mSampDose, resultGUInomScen), ...
                      'matRad:Error');

function test_samplingAnalysisRejectsEvaluationModeMismatch
ct = helper_createAnalysisCt();
cst = helper_createCst();
pln = helper_createAnalysisPlan(ct);
caSampRes = helper_createSamplingResults();
caSampRes(1).evaluationModeBase = 'total';
mSampDose = [1 2 4; ...
             2 3 5];
resultGUInomScen = helper_createNominalResult(ct, cst);

assertExceptionThrown(@() matRad_samplingAnalysis(ct, cst, pln, caSampRes, ...
                                                  mSampDose, resultGUInomScen), ...
                      'matRad:Error');

function test_samplingAnalysisRejectsDoseColumnMismatch
ct = helper_createAnalysisCt();
cst = helper_createCst();
pln = helper_createAnalysisPlan(ct);
caSampRes = helper_createSamplingResults();
mSampDose = [1 2; ...
             2 3];
resultGUInomScen = helper_createNominalResult(ct, cst);

assertExceptionThrown(@() matRad_samplingAnalysis(ct, cst, pln, caSampRes, ...
                                                  mSampDose, resultGUInomScen), ...
                      'matRad:Error');

function test_samplingAnalysisDoesNotReadDeprecatedBioModelQuantities
ct = helper_createAnalysisCt();
cst = helper_createCst();
pln = helper_createAnalysisPlan(ct);
caSampRes = helper_createSamplingResults();
mSampDose = [1 2 4; ...
             2 3 5];
resultGUInomScen = helper_createNominalResult(ct, cst);

analysisText = evalc('matRad_samplingAnalysis(ct, cst, pln, caSampRes, mSampDose, resultGUInomScen);');

assertTrue(isempty(strfind(analysisText, 'Property quantityOpt is deprecated from bioModel')));
assertTrue(isempty(strfind(analysisText, 'Property quantityVis is deprecated from bioModel')));

function test_samplingAnalysisSkipsGammaForPartialSampling
ct = helper_createAnalysisCt();
cst = helper_createCst();
pln = helper_createAnalysisPlan(ct);
caSampRes = helper_createSamplingResults();
mSampDose = [1 2 4; ...
             2 3 5];
resultGUInomScen = helper_createNominalResult(ct, cst);

[~, doseStat] = matRad_samplingAnalysis(ct, cst, pln, caSampRes, ...
                                        mSampDose, resultGUInomScen);

assertEqual(doseStat.gammaAnalysis.status, 'skippedPartialSampling');
assertEqual(doseStat.gammaAnalysis.reason, ...
            'Gamma whole-CT analysis requires sampled statistics for every CT voxel.');
assertTrue(all(isnan(doseStat.gammaAnalysis.gammaCube(:))));
assertTrue(isempty(doseStat.gammaAnalysis.gammaPassRate));
assertEqual(doseStat.expectedDoseDifferenceAnalysis.status, 'computedPartialMask');
assertElementsAlmostEqual(doseStat.expectedDoseDifferenceAnalysis.overReferenceProbabilityCube(1), ...
                          0.8, 'absolute', 1e-12);
assertElementsAlmostEqual(doseStat.expectedDoseDifferenceAnalysis.nearReferenceProbabilityCube(1), ...
                          0.2, 'absolute', 1e-12);
assertTrue(isnan(doseStat.expectedDoseDifferenceAnalysis.overReferenceProbabilityCube(2)));

function test_samplingAnalysisComputesGammaForFullSampling
figureCleaner = onCleanup(@() close('all'));
ct = helper_createFullCoverageCt();
cst = helper_createCst();
pln = helper_createFullCoveragePlan(ct);
caSampRes = helper_createSamplingResults();
mSampDose = ones(prod(ct.cubeDim), numel(caSampRes));
resultGUInomScen = helper_createFullCoverageNominalResult(ct, cst);
existingFig = figure('Visible', 'off');
axes('Parent', existingFig);

[~, doseStat, ~, gammaFig, ~, ~, expectedDoseDifferenceFig] = ...
    matRad_samplingAnalysis(ct, cst, pln, caSampRes, mSampDose, ...
                            resultGUInomScen, 'slice', 1, 'plane', 3, ...
                            'expectedDoseDifferenceDoseWindow', [-2 2]);

assertEqual(doseStat.gammaAnalysis.status, 'computedFullCube');
assertElementsAlmostEqual(doseStat.gammaAnalysis.gammaPassRate, 100, 'absolute', 1e-12);
assertFalse(any(isnan(doseStat.gammaAnalysis.gammaCube(:))));
assertEqual(doseStat.expectedDoseDifferenceAnalysis.status, 'computedFullCube');
assertElementsAlmostEqual(doseStat.expectedDoseDifferenceAnalysis.doseWindow, [-2 2], ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(doseStat.expectedDoseDifferenceAnalysis.nearReferenceProbabilityCube(1), ...
                          1, 'absolute', 1e-12);
assertFalse(isempty(gammaFig));
assertFalse(isequal(gammaFig, existingFig));
assertEqual(get(gammaFig, 'Name'), 'Gamma index analysis');
assertFalse(isempty(expectedDoseDifferenceFig));
helper_closeFigure(gammaFig);
helper_closeFigure(expectedDoseDifferenceFig);

function test_samplingAnalysisStoresRobustnessAnalysisForTargets
ct = helper_createFullCoverageCt();
cst = helper_createTargetCst();
pln = helper_createFullCoveragePlan(ct);
caSampRes = helper_createSamplingResults();
mSampDose = 5 * ones(prod(ct.cubeDim), numel(caSampRes));
resultGUInomScen = helper_createFullCoverageNominalResult(ct, cst);
resultGUInomScen.physicalDose = 5 * ones(ct.cubeDim);

[~, doseStat] = matRad_samplingAnalysis(ct, cst, pln, caSampRes, ...
                                        mSampDose, resultGUInomScen, ...
                                        'robustnessCriteria', [1 1]);

assertTrue(doseStat.robustnessAnalysis.targets(1).isEvaluable);
assertEqual(doseStat.robustnessAnalysis.targets(1).samplingStatus, 'evaluable');
assertElementsAlmostEqual(doseStat.robustnessAnalysis.targets(1).refDose, 5, ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(doseStat.robustnessAnalysis.index1.robPassRate, 100, ...
                          'absolute', 1e-12);

function cst = helper_createCst()
cst = cell(1, 6);
cst{1, 1} = 1;
cst{1, 2} = 'testVoi';
cst{1, 3} = 'OAR';
cst{1, 4} = {[1; 4]};
cst{1, 5} = struct('Visible', true);
cst{1, 6} = {};

function caSampRes = helper_createSamplingResults()
doseGrid = [0 5 10];
volumePoints = [100 50 0; ...
                100 60 0; ...
                100 40 0];
qiValues = [10 20 30];

for i = 1:size(volumePoints, 1)
    caSampRes(i) = helper_createSamplingResult(i, 1, 1, false);
    caSampRes(i).dvh(1).name = 'testVoi';
    caSampRes(i).dvh(1).doseGrid = doseGrid;
    caSampRes(i).dvh(1).volumePoints = volumePoints(i, :);
    caSampRes(i).qi(1).name = 'testVoi';
    caSampRes(i).qi(1).mean = qiValues(i);
end

function caSampRes = helper_createFourDSamplingResults()
caSampRes = helper_createSamplingResults();
ctScenIds = [1 2 1];
for i = 1:3
    caSampRes(i).ctScenId = ctScenIds(i);
    caSampRes(i).scenario.ctScenId = ctScenIds(i);
    caSampRes(i).doseMapping.sourceCtScenId = ctScenIds(i);
    caSampRes(i).doseMapping.referenceCtScenId = 1;
    caSampRes(i).doseMapping.mapped = i == 2;
end

function sample = helper_createSamplingResult(scenarioId, ctScenId, refScen, mapped)
sample = struct();
sample.scenarioId = scenarioId;
sample.ctScenId = ctScenId;
sample.refScen = refScen;
sample.scenario = struct('ctScenId', ctScenId);
sample.doseMapping.sourceCtScenId = ctScenId;
sample.doseMapping.referenceCtScenId = refScen;
sample.doseMapping.method = 'pullDirectDoseMapping';
sample.doseMapping.mapped = mapped;
sample.analysisQuantity = 'physicalDose';
sample.evaluationModeBase = 'perFraction';
sample.dvh = [];
sample.qi = [];

function ct = helper_createCtWithPullDvf()
ct.cubeDim = [2 2 1];
ct.numOfCtScen = 2;
ct.refScen = 1;
ct.resolution.x = 1;
ct.resolution.y = 1;
ct.resolution.z = 1;
ct.dvfMetadata.dvfType = 'pull';
ct.dvfMetadata.dvfUnits = 'voxel';
ct.dvfMetadata.refScen = 1;
ct.dvfMetadata.referenceCtScen = 1;
ct.dvf = cell(1, 2);
ct.dvf{1} = zeros([3 ct.cubeDim]);
ct.dvf{2} = zeros([3 ct.cubeDim]);

function ct = helper_createAnalysisCt()
ct.cubeDim = [2 2 2];
ct.numOfCtScen = 1;
ct.refScen = 1;
ct.resolution.x = 1;
ct.resolution.y = 1;
ct.resolution.z = 1;

function ct = helper_createFourDAnalysisCt()
ct = helper_createAnalysisCt();
ct.numOfCtScen = 2;

function ct = helper_createFullCoverageCt()
ct.cubeDim = [2 2 1];
ct.numOfCtScen = 1;
ct.refScen = 1;
ct.resolution.x = 1;
ct.resolution.y = 1;
ct.resolution.z = 1;

function pln = helper_createAnalysisPlan(ct, ctScenIds)
if nargin < 2
    ctScenIds = ones(3, 1);
end

if isfield(ct, 'numOfCtScen') && any(ctScenIds > ct.numOfCtScen)
    error('Test helper CT scenario ids exceed ct.numOfCtScen.');
end

pln.numOfFractions = 1;
pln.bioModel = matRad_bioModel('photons', 'none');
pln.propOpt.quantityOpt = 'physicalDose';
pln.propOpt.quantityVis = 'physicalDose';
pln.multScen = helper_createSamplingAnalysisModel(ctScenIds);
pln.subIx = [1; 4];

function pln = helper_createFullCoveragePlan(ct)
pln = helper_createAnalysisPlan(ct);
pln.subIx = (1:prod(ct.cubeDim))';

function resultGUI = helper_createNominalResult(ct, cst)
resultGUI.physicalDose = zeros(ct.cubeDim);
resultGUI.physicalDose(1) = 1;
resultGUI.physicalDose(4) = 2;
resultGUI.cst = cst;
resultGUI.analysisQuantity = 'physicalDose';
resultGUI.evaluationModeBase = 'perFraction';

function resultGUI = helper_createFullCoverageNominalResult(ct, cst)
resultGUI.physicalDose = ones(ct.cubeDim);
resultGUI.cst = cst;
resultGUI.analysisQuantity = 'physicalDose';
resultGUI.evaluationModeBase = 'perFraction';

function cst = helper_createTargetCst()
cst = cell(1, 6);
cst{1, 1} = 1;
cst{1, 2} = 'target';
cst{1, 3} = 'TARGET';
cst{1, 4} = {(1:4)'};
cst{1, 5} = struct('Visible', true);
cst{1, 6} = {DoseObjectives.matRad_SquaredDeviation(1, 5)};

function model = helper_createSamplingAnalysisModel(ctScenIds)
ctScenIds = ctScenIds(:);
numScenarios = numel(ctScenIds);
storageSubscripts = [(1:numScenarios)' ones(numScenarios, 1)];
model = matRad_NominalScenario();
components = matRad_createScenarioComponents([1 1 1], 1, 1, ...
                                             {'ct', 'setup', 'range'}, 0);
scenarioValues = zeros(numScenarios, numel(components));
if numScenarios == 3
    scenProb = [0.2; 0.3; 0.5];
else
    scenProb = ones(numScenarios, 1) ./ numScenarios;
end
scenWeight = scenProb;
scenForProb = [ctScenIds scenarioValues];

model.setScenarioRealizations(components, scenarioValues, ctScenIds, ...
                              scenProb, scenWeight, scenForProb, ...
                              storageSubscripts, [numScenarios 1], ...
                              'compact-realization');

function helper_closeFigure(figHandle)
if ~isempty(figHandle) && ishandle(figHandle)
    close(figHandle);
end
