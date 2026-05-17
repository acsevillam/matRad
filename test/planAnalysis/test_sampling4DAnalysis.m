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

sample = helper_createSamplingResult(2, 1, true);

assertEqual(sample.ctScenId, 2);
assertEqual(sample.refScen, 1);
assertTrue(sample.doseMapping.mapped);
assertEqual(sample.doseMapping.sourceCtScenId, 2);
assertEqual(sample.doseMapping.referenceCtScenId, 1);

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
    caSampRes(i) = helper_createSamplingResult(i, 1, false);
    caSampRes(i).dvh(1).name = 'testVoi';
    caSampRes(i).dvh(1).doseGrid = doseGrid;
    caSampRes(i).dvh(1).volumePoints = volumePoints(i, :);
    caSampRes(i).qi(1).name = 'testVoi';
    caSampRes(i).qi(1).mean = qiValues(i);
end

function sample = helper_createSamplingResult(ctScenId, refScen, mapped)
sample = struct();
sample.scenarioId = ctScenId;
sample.ctScenId = ctScenId;
sample.refScen = refScen;
sample.scenario = struct('ctScenId', ctScenId);
sample.doseMapping.sourceCtScenId = ctScenId;
sample.doseMapping.referenceCtScenId = refScen;
sample.doseMapping.method = 'pullDirectDoseMapping';
sample.doseMapping.mapped = mapped;
sample.analysisQuantity = 'physicalDose';
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
ct.resolution.x = 1;
ct.resolution.y = 1;
ct.resolution.z = 1;

function pln = helper_createAnalysisPlan(ct)
pln.numOfFractions = 1;
pln.bioModel = matRad_bioModel('photons', 'none');
pln.propOpt.quantityOpt = 'physicalDose';
pln.propOpt.quantityVis = 'physicalDose';
pln.multScen = helper_createSamplingAnalysisModel();
pln.subIx = [1; 4];

function resultGUI = helper_createNominalResult(ct, cst)
resultGUI.physicalDose = zeros(ct.cubeDim);
resultGUI.physicalDose(1) = 1;
resultGUI.physicalDose(4) = 2;
resultGUI.cst = cst;
resultGUI.analysisQuantity = 'physicalDose';
resultGUI.evaluationModeBase = 'perFraction';

function model = helper_createSamplingAnalysisModel()
storageSubscripts = [(1:3)' ones(3, 1)];
model = matRad_NominalScenario();
components = matRad_createScenarioComponents([1 1 1], 1, 1, ...
                                             {'ct', 'setup', 'range'}, 0);
scenarioValues = zeros(3, numel(components));
ctScenIds = ones(3, 1);
scenProb = [0.2; 0.3; 0.5];
scenWeight = scenProb;
scenForProb = [ctScenIds scenarioValues];

model.setScenarioRealizations(components, scenarioValues, ctScenIds, ...
                              scenProb, scenWeight, scenForProb, ...
                              storageSubscripts, [3 1], ...
                              'compact-realization');
