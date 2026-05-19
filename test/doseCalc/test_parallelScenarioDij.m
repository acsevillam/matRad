function test_suite = test_parallelScenarioDij

test_functions = localfunctions();

initTestSuite;

function test_normalizeScenarioDoseInputConfigDefaultsMemoryLimitToFourGB
ct = struct('refScen', 1);
pln = struct('propOpt', struct());
cfg = ScenarioBatch.Config.matRad_normalizeScenarioDoseInputConfig( ...
    struct(), ct, pln, 'test', MatRad_Config.instance());

assertEqual(cfg.MemoryLimitMB, 4096);

function test_assemblerInsertsMatricesAtOriginalDijIndices
scenarioModel = helper_fixtureScenarioModel();
dij1 = helper_scenarioDij([1 0; 0 1]);
dij2 = helper_scenarioDij([2 0; 0 2]);

dij = ScenarioBatch.Dij.matRad_assembleParallelScenarioDij({dij1; dij2}, [1; 2], scenarioModel);

assertEqual(size(dij.physicalDose), scenarioModel.getDijContainerSize());
assertElementsAlmostEqual(full(dij.physicalDose{1}), ...
                          [1 0; 0 1], 'absolute', 1e-12);
assertTrue(isempty(dij.physicalDose{2}));
assertElementsAlmostEqual(full(dij.physicalDose{3}), ...
                          [2 0; 0 2], 'absolute', 1e-12);
assertEqual(dij.numOfScenarios, 2);
assertEqual(dij.scenarioIds, [1; 2]);
assertEqual(dij.scenarioModel.scenarioIds(), scenarioModel.scenarioIds());
assertEqual(dij.beamNum, dij1.beamNum);

function test_assemblerHandlesAdditionalDoseMatrixFields
scenarioModel = helper_fixtureScenarioModel();
dij1 = helper_scenarioDij([1 0; 0 1]);
dij2 = helper_scenarioDij([2 0; 0 2]);
dij1.mLETDose = {6 * dij1.physicalDose{1}};
dij2.mLETDose = {7 * dij2.physicalDose{1}};
dij1.mAlphaDose = {2 * dij1.physicalDose{1}};
dij2.mAlphaDose = {3 * dij2.physicalDose{1}};
dij1.mSqrtBetaDose = {4 * dij1.physicalDose{1}};
dij2.mSqrtBetaDose = {5 * dij2.physicalDose{1}};
dij1.customQuantity = {8 * dij1.physicalDose{1}};
dij2.customQuantity = {9 * dij2.physicalDose{1}};

dij = ScenarioBatch.Dij.matRad_assembleParallelScenarioDij({dij1; dij2}, [1; 2], scenarioModel);

assertElementsAlmostEqual(full(dij.mLETDose{1}), ...
                          6 * [1 0; 0 1], 'absolute', 1e-12);
assertElementsAlmostEqual(full(dij.mLETDose{3}), ...
                          7 * [2 0; 0 2], 'absolute', 1e-12);
assertElementsAlmostEqual(full(dij.mAlphaDose{1}), ...
                          2 * [1 0; 0 1], 'absolute', 1e-12);
assertElementsAlmostEqual(full(dij.mAlphaDose{3}), ...
                          3 * [2 0; 0 2], 'absolute', 1e-12);
assertElementsAlmostEqual(full(dij.mSqrtBetaDose{1}), ...
                          4 * [1 0; 0 1], 'absolute', 1e-12);
assertElementsAlmostEqual(full(dij.mSqrtBetaDose{3}), ...
                          5 * [2 0; 0 2], 'absolute', 1e-12);
assertElementsAlmostEqual(full(dij.customQuantity{1}), ...
                          8 * [1 0; 0 1], 'absolute', 1e-12);
assertElementsAlmostEqual(full(dij.customQuantity{3}), ...
                          9 * [2 0; 0 2], 'absolute', 1e-12);

function test_assemblerInsertsBiologicalCtMetadata
scenarioModel = helper_fixtureTwoCtScenarioModel();
dij1 = helper_scenarioDij([1 0; 0 1]);
dij2 = helper_scenarioDij([2 0; 0 2]);
dij1 = helper_addBiologicalCtMetadata(dij1, [11; 12], [], [21; 22], [], ...
                                      [31; 32], [], [41; 42], [], [1; 2], [], [3; 4], []);
dij2 = helper_addBiologicalCtMetadata(dij2, [], [13; 14], [], [23; 24], ...
                                      [], [33; 34], [], [43; 44], [], [1; 2], [], [7; 8]);

dij = ScenarioBatch.Dij.matRad_assembleParallelScenarioDij({dij1; dij2}, [1; 2], ...
                                                           scenarioModel);

assertEqual(dij.ax{1}, [11; 12]);
assertEqual(dij.ax{2}, [13; 14]);
assertEqual(dij.bx{1}, [21; 22]);
assertEqual(dij.bx{2}, [23; 24]);
assertEqual(dij.abx{1}, [31; 32]);
assertEqual(dij.abx{2}, [33; 34]);
assertEqual(dij.gamma{1}, [41; 42]);
assertEqual(dij.gamma{2}, [43; 44]);
assertEqual(dij.ixDose{1}, [1; 2]);
assertEqual(dij.ixDose{2}, [1; 2]);
assertEqual(dij.vTissueIndex{1}, [3; 4]);
assertEqual(dij.vTissueIndex{2}, [7; 8]);

function test_assemblerAcceptsIdenticalBiologicalMetadataForSharedCt
scenarioModel = helper_fixtureScenarioModel();
dij1 = helper_scenarioDij([1 0; 0 1]);
dij2 = helper_scenarioDij([2 0; 0 2]);
dij1.ax = {[1; 2]};
dij2.ax = {[1; 2]};
dij1.vTissueIndex = {[3; 4]};
dij2.vTissueIndex = {[3; 4]};

dij = ScenarioBatch.Dij.matRad_assembleParallelScenarioDij({dij1; dij2}, [1; 2], ...
                                                           scenarioModel);

assertEqual(dij.ax{1}, [1; 2]);
assertEqual(dij.vTissueIndex{1}, [3; 4]);

function test_assemblerRejectsDivergentBiologicalMetadataForSharedCt
scenarioModel = helper_fixtureScenarioModel();
dij1 = helper_scenarioDij([1 0; 0 1]);
dij2 = helper_scenarioDij([2 0; 0 2]);
dij1.ax = {[1; 2]};
dij2.ax = {[2; 3]};

assertExceptionThrown(@() ScenarioBatch.Dij.matRad_assembleParallelScenarioDij( ...
                                                                               {dij1; dij2}, [1; 2], scenarioModel), 'matRad:Error');

function test_assemblerRejectsBiologicalMetadataWithWrongLength
scenarioModel = helper_fixtureScenarioModel();
dij1 = helper_scenarioDij([1 0; 0 1]);
dij2 = helper_scenarioDij([2 0; 0 2]);
dij1.ax = {[1; 2]};
dij2.ax = {[1; 2; 3]};

assertExceptionThrown(@() ScenarioBatch.Dij.matRad_assembleParallelScenarioDij( ...
                                                                               {dij1; dij2}, [1; 2], scenarioModel), 'matRad:Error');

function test_assemblerRejectsBiologicalMetadataWithInvalidIndices
scenarioModel = helper_fixtureScenarioModel();
dij1 = helper_scenarioDij([1 0; 0 1]);
dij2 = helper_scenarioDij([2 0; 0 2]);
dij1.ixDose = {[1; 2]};
dij2.ixDose = {[1; 3]};

assertExceptionThrown(@() ScenarioBatch.Dij.matRad_assembleParallelScenarioDij( ...
                                                                               {dij1; dij2}, [1; 2], scenarioModel), 'matRad:Error');

function test_assemblerRejectsIncompatibleGrids
scenarioModel = helper_fixtureScenarioModel();
dij1 = helper_scenarioDij([1 0; 0 1]);
dij2 = helper_scenarioDij([2 0; 0 2]);
dij2.doseGrid.numOfVoxels = 3;

assertExceptionThrown(@() ScenarioBatch.Dij.matRad_assembleParallelScenarioDij( ...
                                                                               {dij1; dij2}, [1; 2], scenarioModel), 'matRad:Error');

function test_assemblerRejectsIncompatibleBixelCount
scenarioModel = helper_fixtureScenarioModel();
dij1 = helper_scenarioDij([1 0; 0 1]);
dij2 = helper_scenarioDij([2 0; 0 2]);
dij2.totalNumOfBixels = 3;

assertExceptionThrown(@() ScenarioBatch.Dij.matRad_assembleParallelScenarioDij( ...
                                                                               {dij1; dij2}, [1; 2], scenarioModel), 'matRad:Error');

function test_assemblerRejectsDuplicateScenarioIds
scenarioModel = helper_fixtureScenarioModel();
dij1 = helper_scenarioDij([1 0; 0 1]);
dij2 = helper_scenarioDij([2 0; 0 2]);

assertExceptionThrown(@() ScenarioBatch.Dij.matRad_assembleParallelScenarioDij( ...
                                                                               {dij1; dij2}, [1; 1], scenarioModel), 'matRad:Error');

function test_assemblerRejectsIncompleteScenarioIds
scenarioModel = helper_fixtureScenarioModel();
dij1 = helper_scenarioDij([1 0; 0 1]);

assertExceptionThrown(@() ScenarioBatch.Dij.matRad_assembleParallelScenarioDij( ...
                                                                               {dij1}, 1, scenarioModel), 'matRad:Error');

function test_parallelPencilBeamMultiscenMatchesSerialOrFallback
[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);

pln.propDoseCalc.UseParallel = false;
dijSerial = matRad_calcDoseInfluence(ct, cst, stf, pln);

pln.propDoseCalc.UseParallel = true;
dijParallel = matRad_calcDoseInfluence(ct, cst, stf, pln);

helper_assertDijAlmostEqual(dijSerial, dijParallel);

function test_parallelPencilBeamMultiscenUsesParforWhenAvailable
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is unavailable.');
end

[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);
pln.propDoseCalc.enableDijSampling = false;

pln.propDoseCalc.UseParallel = false;
dijSerial = matRad_calcDoseInfluence(ct, cst, stf, pln);

pln.propDoseCalc.UseParallel = true;
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
[dijParallel, useParallel] = ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine);

if ~useParallel
    moxunit_throw_test_skipped_exception(['Parallel scenario dij ', ...
                                          'calculation could not be activated in this environment.']);
end
helper_assertDijAlmostEqual(dijSerial, dijParallel);

function test_parallelPencilBeamMultiscenRespectsWorkerUpperBound
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is unavailable.');
end

[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);
pln.propDoseCalc.UseParallel = true;
pln.propDoseCalc.enableDijSampling = false;
pln.propDoseCalc.parallelOptions = struct('workerUpperBound', 1);
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

[dij, useParallel] = ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine);

assertFalse(useParallel);
assertTrue(isempty(dij));

function test_parallelPencilBeamMultiscenIncreasesSmallPoolWhenSafe
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception(['Parallel Computing ', ...
                                          'Toolbox is unavailable.']);
end
if helper_maxParallelPoolWorkers() < 2
    moxunit_throw_test_skipped_exception(['Local parallel pool ', ...
                                          'supports fewer than two workers.']);
end

cleanup = helper_preserveParallelPool(); %#ok<NASGU>
matRad_configureParallelPoolSize(1, 'parallel scenario dij test', ...
                                 MatRad_Config.instance());

[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);
pln.propDoseCalc.UseParallel = true;
pln.propDoseCalc.enableDijSampling = false;
pln.propDoseCalc.parallelOptions = struct('workerUpperBound', 2);
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

[dij, useParallel] = ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine);

if ~useParallel
    moxunit_throw_test_skipped_exception(['Parallel scenario dij ', ...
                                          'calculation could not be activated in this environment.']);
end
pPool = gcp('nocreate');
assertTrue(~isempty(pPool));
assertTrue(pPool.NumWorkers > 1);
assertTrue(pPool.NumWorkers <= 2);
assertFalse(isempty(dij));

function test_parallelPencilBeamMultiscenRejectsUnknownParallelOption
[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);
pln.propDoseCalc.UseParallel = true;
pln.propDoseCalc.enableDijSampling = false;
pln.propDoseCalc.parallelOptions = struct('workers', 2);
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

assertExceptionThrown(@() ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine), 'matRad:Error');

function test_parallelPencilBeamMultiscenRejectsFractionalWorkerUpperBound
[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);
pln.propDoseCalc.UseParallel = true;
pln.propDoseCalc.enableDijSampling = false;
pln.propDoseCalc.parallelOptions = struct('workerUpperBound', 1.5);
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

assertExceptionThrown(@() ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine), 'matRad:Error');

function test_parallelPencilBeamMultiscenEngineHandleIsWorkerSafe
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is unavailable.');
end

[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);
pln.propDoseCalc.enableDijSampling = false;

pln.propDoseCalc.UseParallel = false;
dijSerial = matRad_calcDoseInfluence(ct, cst, stf, pln);

engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
engine.UseParallel = true;
engine.randomSeed = 19;
originalRandomSeed = engine.randomSeed;
originalScenarioIds = engine.multScen.scenarioIds();

pln.propDoseCalc = engine;
[dijParallel, useParallel] = ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine);

assertTrue(engine.UseParallel);
assertEqual(engine.randomSeed, originalRandomSeed);
assertEqual(engine.multScen.scenarioIds(), originalScenarioIds);
if ~useParallel
    moxunit_throw_test_skipped_exception(['Parallel scenario dij ', ...
                                          'calculation could not be activated in this environment.']);
end
helper_assertDijAlmostEqual(dijSerial, dijParallel);

function test_parallelPencilBeamMultiscenFallsBackForStochasticSampling
[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);
pln.propDoseCalc.UseParallel = true;
pln.propDoseCalc.enableDijSampling = true;
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

[dij, useParallel] = ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine);

assertFalse(useParallel);
assertTrue(isempty(dij));

function test_parallelSupportQueryRejectsStochasticPhotonDij
[~, ~, pln, ~] = helper_photonTestDataFixture();
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

[isSupported, reason] = ScenarioBatch.Parallel.matRad_supportsParallelScenarioDij(engine);

assertFalse(isSupported);
assertEqual(reason, 'stochasticDijSampling');

function test_parallelSupportQueryAcceptsAnalyticalParticleDij
[~, ~, pln, ~] = helper_particleBioTestDataFixture();
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

[isSupported, reason] = ScenarioBatch.Parallel.matRad_supportsParallelScenarioDij(engine);

assertTrue(isSupported);
assertEqual(reason, '');

function test_parallelParticleBioMultiscenMatchesSerialWhenAvailable
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is unavailable.');
end

[ct, cst, pln, stf] = helper_particleBioTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);

pln.propDoseCalc.UseParallel = false;
dijSerial = matRad_calcDoseInfluence(ct, cst, stf, pln);

pln.propDoseCalc.UseParallel = true;
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
[dijParallel, useParallel] = ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine);

if ~useParallel
    moxunit_throw_test_skipped_exception(['Parallel particle dij ', ...
                                          'calculation could not be activated in this environment.']);
end
helper_assertDijAlmostEqual(dijSerial, dijParallel);

function test_parallelParticleLemMultiscenMatchesSerialWhenAvailable
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is unavailable.');
end

[ct, cst, pln, stf] = helper_particleLemTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);

pln.propDoseCalc.UseParallel = false;
dijSerial = matRad_calcDoseInfluence(ct, cst, stf, pln);

pln.propDoseCalc.UseParallel = true;
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
[dijParallel, useParallel] = ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine);

if ~useParallel
    moxunit_throw_test_skipped_exception(['Parallel particle LEM dij ', ...
                                          'calculation could not be activated in this environment.']);
end
helper_assertDijAlmostEqual(dijSerial, dijParallel);

function test_parallelHelperFallsBackForSingleScenario
[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = matRad_NominalScenario(ct);
pln.propDoseCalc.UseParallel = true;
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

[dij, useParallel] = ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine);

assertFalse(useParallel);
assertTrue(isempty(dij));

function test_parallelHelperFallsBackForNonPencilBeamEngine
pln = struct('radiationMode', 'brachy', 'machine', 'Generic');
pln.propDoseCalc.UseParallel = true;
engine = DoseEngines.matRad_TG43BrachyEngine(pln);

[dij, useParallel] = ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios([], [], [], pln, engine);

assertFalse(useParallel);
assertTrue(isempty(dij));

function test_parallelHelperFallsBackWithoutParallelToolbox
if helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is available.');
end

[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);
pln.propDoseCalc.UseParallel = true;
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

[dij, useParallel] = ScenarioBatch.Parallel.matRad_calcDoseInfluenceParallelScenarios(ct, cst, stf, pln, engine);

assertFalse(useParallel);
assertTrue(isempty(dij));

function scenarioModel = helper_fixtureScenarioModel()
ct.numOfCtScen = 1;
ct.cubeDim = [2 1 1];
ct.resolution = struct('x', 1, 'y', 1, 'z', 1);
scenarioModel = matRad_NominalScenario(ct);
dimensions = matRad_createScenarioComponents([1 1 1], 1, 1);
scenarioValues = [0 0 0 0 0; 1 0 0 0 0];
ctScenIds = [1; 1];
storageSubscripts = [1 1 1; 1 3 1];
storageSize = [1 3 1];
scenarioWeights = [0.5; 0.5];
scenForProb = [ctScenIds scenarioValues];
scenarioModel.setScenarioRealizations(dimensions, scenarioValues, ctScenIds, ...
                                      scenarioWeights, scenarioWeights, scenForProb, storageSubscripts, storageSize, ...
                                      'legacy-grid');

function scenarioModel = helper_fixtureTwoCtScenarioModel()
ct.numOfCtScen = 2;
ct.cubeDim = [2 1 1];
ct.resolution = struct('x', 1, 'y', 1, 'z', 1);
scenarioModel = matRad_NominalScenario(ct);
dimensions = matRad_createScenarioComponents([1 1 1], 1, 1);
scenarioValues = [0 0 0 0 0; 0 0 0 0 0];
ctScenIds = [1; 2];
storageSubscripts = [1 1 1; 2 1 1];
storageSize = [2 1 1];
scenarioWeights = [0.5; 0.5];
scenForProb = [ctScenIds scenarioValues];
scenarioModel.setScenarioRealizations(dimensions, scenarioValues, ctScenIds, ...
                                      scenarioWeights, scenarioWeights, scenForProb, storageSubscripts, storageSize, ...
                                      'legacy-grid');

function dij = helper_scenarioDij(matrix)
matrix = sparse(matrix);
dij.ctGrid = helper_doseGrid(size(matrix, 1));
dij.doseGrid = dij.ctGrid;
dij.numOfBeams = 1;
dij.numOfScenarios = 1;
dij.scenarioModel = matRad_NominalScenario();
dij.scenarioIds = 1;
dij.numOfRaysPerBeam = size(matrix, 2);
dij.totalNumOfBixels = size(matrix, 2);
dij.totalNumOfRays = size(matrix, 2);
dij.beamNum = ones(size(matrix, 2), 1);
dij.rayNum = (1:size(matrix, 2))';
dij.bixelNum = ones(size(matrix, 2), 1);
dij.minMU = zeros(size(matrix, 2), 1);
dij.maxMU = inf(size(matrix, 2), 1);
dij.numOfParticlesPerMU = 1e6 * ones(size(matrix, 2), 1);
dij.physicalDose = {matrix};

function grid = helper_doseGrid(numVoxels)
grid.dimensions = [numVoxels 1 1];
grid.numOfVoxels = numVoxels;
grid.resolution = struct('x', 1, 'y', 1, 'z', 1);
grid.x = 1;
grid.y = 1:numVoxels;
grid.z = 1;

function [ct, cst, pln, stf] = helper_photonTestDataFixture()
testDataPath = fullfile(fileparts(mfilename('fullpath')), ...
                        '..', 'testData', 'photons_testData.mat');
data = load(testDataPath, 'ct', 'cst', 'pln', 'stf');
ct = data.ct;
cst = data.cst;
pln = data.pln;
stf = data.stf;
pln.propDoseCalc.engine = 'SVDPB';

function [ct, cst, pln, stf] = helper_particleBioTestDataFixture()
testDataPath = fullfile(fileparts(mfilename('fullpath')), ...
                        '..', 'testData', 'protons_testData.mat');
data = load(testDataPath, 'ct', 'cst', 'pln', 'stf');
ct = data.ct;
cst = data.cst;
pln = data.pln;
stf = data.stf;
pln.propDoseCalc.engine = 'HongPB';
pln.bioParam = matRad_bioModel(pln.radiationMode, 'RBExD', 'MCN');

function [ct, cst, pln, stf] = helper_particleLemTestDataFixture()
testDataPath = fullfile(fileparts(mfilename('fullpath')), ...
                        '..', 'testData', 'carbon_testData.mat');
data = load(testDataPath, 'ct', 'cst', 'pln', 'stf');
ct = data.ct;
cst = data.cst;
pln = data.pln;
stf = data.stf;
pln.propDoseCalc.engine = 'HongPB';
pln.bioParam = matRad_bioModel(pln.radiationMode, 'RBExD', 'LEM');

function dij = helper_addBiologicalCtMetadata(dij, ax1, ax2, bx1, bx2, abx1, abx2, ...
                                              gamma1, gamma2, ixDose1, ixDose2, vTissueIndex1, vTissueIndex2)
dij.ax = {ax1; ax2};
dij.bx = {bx1; bx2};
dij.abx = {abx1; abx2};
dij.gamma = {gamma1; gamma2};
dij.ixDose = {ixDose1; ixDose2};
dij.vTissueIndex = {vTissueIndex1; vTissueIndex2};

function scenario = helper_rangeRandomScenario(ct)
scenario = matRad_RandomScenarios(ct);
scenario.nSamples = 2;
scenario.randomSeed = 31;
scenario.scenarioDimensionActive = {'ct', 'range'};

function helper_assertDijAlmostEqual(expectedDij, actualDij)
assertEqual(actualDij.scenarioIds, expectedDij.scenarioIds);
assertEqual(actualDij.numOfScenarios, expectedDij.numOfScenarios);
assertEqual(actualDij.totalNumOfBixels, expectedDij.totalNumOfBixels);
assertEqual(actualDij.beamNum, expectedDij.beamNum);
assertEqual(actualDij.rayNum, expectedDij.rayNum);
assertEqual(actualDij.bixelNum, expectedDij.bixelNum);
matrixFields = {'physicalDose', 'mLETDose', 'mAlphaDose', 'mSqrtBetaDose'};
for fieldIx = 1:numel(matrixFields)
    helper_assertDijCellFieldAlmostEqual(expectedDij, actualDij, ...
                                         matrixFields{fieldIx}, 1e-10);
end

bioFields = {'ax', 'bx', 'abx', 'gamma', 'ixDose', 'vTissueIndex'};
for fieldIx = 1:numel(bioFields)
    helper_assertDijCellFieldAlmostEqual(expectedDij, actualDij, ...
                                         bioFields{fieldIx}, 1e-12);
end

function helper_assertDijCellFieldAlmostEqual(expectedDij, actualDij, fieldName, tolerance)
if ~isfield(expectedDij, fieldName)
    assertFalse(isfield(actualDij, fieldName));
    return
end

assertTrue(isfield(actualDij, fieldName));
assertEqual(size(actualDij.(fieldName)), size(expectedDij.(fieldName)));

for i = 1:numel(expectedDij.(fieldName))
    if isempty(expectedDij.(fieldName){i})
        assertTrue(isempty(actualDij.(fieldName){i}));
    else
        assertElementsAlmostEqual(full(actualDij.(fieldName){i}), ...
                                  full(expectedDij.(fieldName){i}), 'absolute', tolerance);
    end
end

function available = helper_parallelComputingAvailable()
available = false;
if exist('parpool', 'file') ~= 2 || exist('gcp', 'file') ~= 2
    return
end

try
    [available, ~] = license('checkout', 'Distrib_Computing_Toolbox');
catch
    available = false;
end

if isempty(available)
    available = false;
end
available = logical(available);

function maxWorkers = helper_maxParallelPoolWorkers()
maxWorkers = 0;
try
    cluster = parcluster();
    maxWorkers = cluster.NumWorkers;
catch
end

if isempty(maxWorkers) || ~isnumeric(maxWorkers) || maxWorkers < 1
    try
        maxWorkers = feature('numcores');
    catch
        maxWorkers = 0;
    end
end

function cleanup = helper_preserveParallelPool()
try
    pPool = gcp('nocreate');
catch
    pPool = [];
end
if isempty(pPool)
    originalNumWorkers = [];
else
    originalNumWorkers = pPool.NumWorkers;
end
cleanup = onCleanup(@() helper_restoreParallelPool(originalNumWorkers));

function helper_restoreParallelPool(originalNumWorkers)
try
    pPool = gcp('nocreate');
catch
    pPool = [];
end
if ~isempty(pPool)
    delete(pPool);
end
if ~isempty(originalNumWorkers)
    try
        parpool(originalNumWorkers);
    catch
    end
end
