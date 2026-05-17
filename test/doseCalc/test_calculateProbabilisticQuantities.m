function test_suite = test_calculateProbabilisticQuantities

test_functions = localfunctions();

initTestSuite;

function test_probPhysicalDoseExpectedAndOmega
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();

[plnOut, dijProbContext] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);
dijProb = plnOut.propOpt.dij_prob;

expected = [2.5 0; 0 3.5; 1.75 1; 0 4];
assertElementsAlmostEqual(full(dijProb.expected), expected, 'absolute', 1e-12);
assertElementsAlmostEqual(full(dijProb.Omega{1}), diag([0.75 0.75]), ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(full(dijProb.Omega{2}), ...
                          [0.1875 0; 0 3], 'absolute', 1e-12);
assertEqual(dijProb.quantity, 'physicalDose');
assertEqual(dijProb.quantityField, 'physicalDose');
assertEqual(dijProb.probabilisticMode, 'PROB');
assertElementsAlmostEqual(dijProb.scenarioWeights, [0.25; 0.75], ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(full(dijProbContext.physicalDose{1}), ...
                          full(dij.physicalDose{1}), 'absolute', 1e-12);
assertTrue(any(abs(full(dijProbContext.physicalDose{1}(:)) - ...
                   expected(:)) > 1e-12));
assertEqual(dijProbContext.numOfScenarios, 1);
assertTrue(isa(dijProbContext.scenarioModel, 'matRad_NominalScenario'));
assertEqual(dijProbContext.probabilisticQuantity, 'physicalDose');
assertEqual(dijProbContext.probabilisticQuantityField, 'physicalDose');

function test_probBatchSizeDoesNotChangeResult
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.BatchSize = 1;
[plnBatch, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

cfg.BatchSize = 99;
[plnFull, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

assertElementsAlmostEqual(full(plnBatch.propOpt.dij_prob.expected), ...
                          full(plnFull.propOpt.dij_prob.expected), 'absolute', 1e-12);
assertElementsAlmostEqual(full(plnBatch.propOpt.dij_prob.Omega{2}), ...
                          full(plnFull.propOpt.dij_prob.Omega{2}), 'absolute', 1e-12);

function test_probResponseIsConsistentWithSelectedVois
[ct, cst, pln, dij, cfg] = helper_partialSelectionFixture();

[plnOut, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);
dijProb = plnOut.propOpt.dij_prob;

expected = [2.5 0; 0 3.5; 1.75 1; 0 0];
assertElementsAlmostEqual(full(dijProb.expected), expected, 'absolute', 1e-12);
assertEqual(dijProb.voiSubIx{1}, [1; 2]);
assertEqual(dijProb.voiSubIx{2}, 3);
assertTrue(isempty(dijProb.voiSubIx{3}));
assertElementsAlmostEqual(full(dijProb.Omega{1}), diag([0.75 0.75]), ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(full(dijProb.Omega{2}), ...
                          [0.1875 0; 0 0], 'absolute', 1e-12);
assertTrue(isempty(dijProb.Omega{3}));

function test_probVoiRowsFollowOverlapPriorities
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cst{1, 5}.Priority = 1;
cst{2, 5}.Priority = 2;
cst{1, 6} = {struct()};
cst{2, 4}{1} = [2; 3; 4];
cst{2, 6} = {struct()};

[plnOut, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);
dijProb = plnOut.propOpt.dij_prob;

assertEqual(dijProb.voiSubIx{1}, [1; 2]);
assertEqual(dijProb.voiSubIx{2}, [3; 4]);

function test_probConstRbeScalesExpectedAndOmega
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
dij.RBE = 1.1;
cfg.Quantity = 'RBExD';

[plnOut, dijProbContext] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);
dijProb = plnOut.propOpt.dij_prob;

assertElementsAlmostEqual(full(dijProb.expected(1:2, :)), ...
                          1.1 * [2.5 0; 0 3.5], 'absolute', 1e-12);
assertElementsAlmostEqual(full(dijProb.Omega{1}), ...
                          1.21 * diag([0.75 0.75]), 'absolute', 1e-12);
assertEqual(dijProb.quantity, 'RBExD');
assertEqual(dijProbContext.RBE, 1);

function test_probAcceptsExplicitLinearQuantityField
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
dij.effectDose = cell(size(dij.physicalDose));
dij.effectDose{1} = 2 .* dij.physicalDose{1};
dij.effectDose{2} = 2 .* dij.physicalDose{2};
cfg.Quantity = 'effect';
cfg.QuantityField = 'effectDose';

[plnOut, dijProbContext] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);
dijProb = plnOut.propOpt.dij_prob;

assertElementsAlmostEqual(full(dijProb.expected(1:2, :)), ...
                          2 .* [2.5 0; 0 3.5], 'absolute', 1e-12);
assertElementsAlmostEqual(full(dijProb.Omega{1}), ...
                          4 .* diag([0.75 0.75]), 'absolute', 1e-12);
assertEqual(dijProb.quantity, 'effect');
assertEqual(dijProb.quantityField, 'effectDose');
assertEqual(dijProbContext.probabilisticQuantity, 'effect');
assertEqual(dijProbContext.probabilisticQuantityField, 'effectDose');

function test_probOmegaMatchesCenteredScenarioAccumulation
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
dij.totalNumOfBixels = 3;
dij.physicalDose{1} = sparse([1 0 5; 0 2 7; 1 1 3; 0 1 9]);
dij.physicalDose{2} = sparse([3 0 5; 0 4 7; 2 1 3; 0 5 9]);

[plnOut, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);
dijProb = plnOut.propOpt.dij_prob;
targetRows = cst{1, 4}{1};
expectedRows = dijProb.expected(targetRows, :);
scenarioRows = {dij.physicalDose{1}(targetRows, :), ...
                dij.physicalDose{2}(targetRows, :)};
expectedOmega = helper_manualCenteredOmega(scenarioRows, ...
                                           dijProb.scenarioWeights, expectedRows, dij.totalNumOfBixels);

assertEqual(size(dijProb.Omega{1}), [3 3]);
assertElementsAlmostEqual(full(dijProb.Omega{1}), ...
                          full(expectedOmega), 'absolute', 1e-12);
assertEqual(nnz(dijProb.Omega{1}(:, 3)), 0);
assertEqual(nnz(dijProb.Omega{1}(3, :)), 0);

function test_probScenarioBatchDiskMatchesRecompute
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
[plnDisk, dijDiskContext] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

cfg.SecondPassStrategy = 'recompute';
[plnRecompute, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

helper_assertProbAlmostEqual(plnDisk.propOpt.dij_prob, ...
                             plnRecompute.propOpt.dij_prob);
assertEqual(plnDisk.propOpt.dij_prob.precomputeMode, 'scenario-batch');
assertEqual(plnDisk.propOpt.dij_prob.secondPassStrategy, 'disk');
assertEqual(plnRecompute.propOpt.dij_prob.secondPassStrategy, 'recompute');
assertEqual(dijDiskContext.numOfScenarios, 1);
assertElementsAlmostEqual(full(dijDiskContext.physicalDose{1}), ...
                          full(dij.physicalDose{1}), 'absolute', 1e-12);
assertTrue(any(abs(full(dijDiskContext.physicalDose{1}(:)) - ...
                   full(plnDisk.propOpt.dij_prob.expected(:))) > 1e-12));

function test_probScenarioBatchPartialSelectionDiskMatchesRecompute
[ct, cst, pln, dij, cfg] = helper_partialSelectionFixture();
cfg.BatchSize = 1;

[plnDisk, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

assertEqual(plnDisk.propOpt.dij_prob.secondPassStrategy, 'disk');

cfg.SecondPassStrategy = 'recompute';
[plnRecompute, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

helper_assertProbAlmostEqual(plnDisk.propOpt.dij_prob, ...
                             plnRecompute.propOpt.dij_prob);
assertEqual(plnRecompute.propOpt.dij_prob.secondPassStrategy, ...
            'recompute');

function test_probScenarioBatchSerialDetailedProgressIsAggregated
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.BatchSize = 1;
cfg.SecondPassStrategy = 'recompute';
cfg.ProgressLevel = 'detailed';

[logCleanup, startLogCount] = helper_captureMatRadLog();
matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

messages = helper_matRadLogMessages(startLogCount);
assertTrue(helper_anyMessageContains(messages, 'Expected influence progress 8/8'));
assertTrue(helper_anyMessageContains(messages, 'scenario 2/2, batch 4/4'));
assertTrue(helper_anyMessageContains(messages, 'omega progress 8/8'));
assertTrue(helper_anyMessageContains(messages, 'scenario 2/2, VOI 2, batch 2/2'));

function test_probScenarioBatchUseParallelPrecomputedDiskMatchesSerial
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.BatchSize = 1;

cfg.UseParallel = false;
[plnSerial, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

cfg.UseParallel = true;
[plnParallel, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

helper_assertProbAlmostEqual(plnSerial.propOpt.dij_prob, ...
                             plnParallel.propOpt.dij_prob);
helper_assertScenarioBatchSizeDisk(plnParallel.propOpt.dij_prob);

function test_probScenarioBatchUseParallelPrecomputedRecomputeMatchesSerial
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.BatchSize = 1;
cfg.SecondPassStrategy = 'recompute';

cfg.UseParallel = false;
[plnSerial, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

cfg.UseParallel = true;
[plnParallel, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

helper_assertProbAlmostEqual(plnSerial.propOpt.dij_prob, ...
                             plnParallel.propOpt.dij_prob);
helper_assertScenarioBatchSizeRecompute(plnParallel.propOpt.dij_prob);

function test_probScenarioBatchAcceptsDijWithoutCfg
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();

[plnDefault, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij);
[plnConfigured, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

helper_assertProbAlmostEqual(plnDefault.propOpt.dij_prob, ...
                             plnConfigured.propOpt.dij_prob);

function test_probScenarioBatchRejectsDuplicatePrecomputedDijInputs
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.PrecomputedDij = dij;

assertExceptionThrown(@() matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg), ...
                      'matRad:Error');

function test_probScenarioBatchRecomputesScenarioDijWithoutDijArgument
[ct, cst, pln, stf] = helper_photonTestDataFixture();
cfg.SecondPassStrategy = 'disk';
cfg.KeepCache = false;
cfg.BatchSize = 10000;
cfg.targetStructSel = 1;

[plnResult, dijContext] = matRad_calculateProbabilisticQuantities(ct, cst, stf, pln, cfg);
dijProb = plnResult.propOpt.dij_prob;

assertEqual(size(dijProb.expected, 1), dijContext.doseGrid.numOfVoxels);
assertEqual(size(dijProb.expected, 2), dijContext.totalNumOfBixels);
assertTrue(nnz(dijProb.expected) > 0);
assertFalse(isfield(dijProb, 'cacheDir'));
helper_assertScenarioBatchSizeDisk(dijProb);
assertEqual(dijContext.numOfScenarios, 1);

function test_probScenarioBatchRecomputeDoesNotMutateEngineHandle
[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = matRad_NominalScenario(ct);
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
engine.UseParallel = true;
engine.doseGrid.resolution = struct('x', 5, 'y', 6, 'z', 7);
engine.selectVoxelsInScenarios = 'all';
engine.randomSeed = 17;
originalScenarioIds = engine.multScen.scenarioIds();
originalRandomSeed = engine.randomSeed;
originalDoseGrid = engine.doseGrid;
originalSelectVoxelsInScenarios = engine.selectVoxelsInScenarios;
pln.propDoseCalc = engine;

cfg.SecondPassStrategy = 'recompute';
cfg.KeepCache = false;
cfg.BatchSize = 10000;
cfg.targetStructSel = 1;

[~, dijContext] = matRad_calculateProbabilisticQuantities(ct, cst, stf, ...
                                                          pln, cfg);

assertTrue(engine.UseParallel);
assertEqual(engine.multScen.scenarioIds(), originalScenarioIds);
assertEqual(engine.randomSeed, originalRandomSeed);
assertEqual(engine.selectVoxelsInScenarios, ...
            originalSelectVoxelsInScenarios);
assertEqual(engine.doseGrid, originalDoseGrid);
assertEqual(dijContext.doseGrid.resolution, ...
            engine.doseGrid.resolution);

function test_probScenarioBatchCollectTimingReportsRealParallelWhenAvailable
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is unavailable.');
end

[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
engine.UseParallel = true;
engine.randomSeed = 29;
originalRandomSeed = engine.randomSeed;
originalScenarioIds = engine.multScen.scenarioIds();
pln.propDoseCalc = engine;

cfg.SecondPassStrategy = 'recompute';
cfg.KeepCache = false;
cfg.BatchSize = 10000;
cfg.targetStructSel = 1;
cfg.UseParallel = true;
cfg.CollectTiming = true;
cfg.ProgressLevel = 'detailed';

[logCleanup, startLogCount] = helper_captureMatRadLog();
[plnResult, ~] = matRad_calculateProbabilisticQuantities(ct, cst, stf, pln, cfg);
timing = plnResult.propOpt.dij_prob.timing;

assertTrue(engine.UseParallel);
assertEqual(engine.randomSeed, originalRandomSeed);
assertEqual(engine.multScen.scenarioIds(), originalScenarioIds);
assertEqual(timing.calculationMode, 'PROB');
if ~timing.parallelScenario.firstPass || ~timing.parallelScenario.omega
    moxunit_throw_test_skipped_exception(['Parallel probabilistic scenario-batch ', ...
                                          'could not be activated in this environment.']);
end
assertTrue(timing.parallelScenario.firstPass);
assertTrue(timing.parallelScenario.omega);
messages = helper_matRadLogMessages(startLogCount);
assertTrue(helper_anyMessageContains(messages, 'Expected influence progress'));
assertTrue(helper_anyMessageContains(messages, 'omega progress'));

function test_probScenarioBatchSmallBatchesDiskMatchRecompute
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.BatchSize = 1;
[plnDisk, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

cfg.SecondPassStrategy = 'recompute';
[plnResult, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

helper_assertProbAlmostEqual(plnDisk.propOpt.dij_prob, ...
                             plnResult.propOpt.dij_prob);
assertEqual(plnResult.propOpt.dij_prob.secondPassStrategy, 'recompute');
helper_assertScenarioBatchSizeRecompute(plnResult.propOpt.dij_prob);

function test_probScenarioBatchDiskCacheKeepsDistinctHashFolders
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cacheRoot = tempname();
cleanup = onCleanup(@() helper_deleteDirIfExists(cacheRoot));
cfg.CacheRoot = cacheRoot;
cfg.KeepCache = true;

[plnRun1, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);
[plnRun2, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

prob_1 = plnRun1.propOpt.dij_prob;
prob_2 = plnRun2.propOpt.dij_prob;
assertTrue(isfield(prob_1, 'cacheDir'));
assertTrue(isfield(prob_2, 'cacheDir'));
helper_assertScenarioBatchSizeDisk(prob_1);
helper_assertScenarioBatchSizeDisk(prob_2);
assertFalse(strcmp(prob_1.cacheDir, prob_2.cacheDir));
assertEqual(exist(prob_1.cacheDir, 'dir'), 7);
assertEqual(exist(prob_2.cacheDir, 'dir'), 7);
assertEqual(exist(fullfile(prob_1.cacheDir, 'metadata.mat'), 'file'), 2);
metadataData = load(fullfile(prob_1.cacheDir, 'metadata.mat'));
assertTrue(isfield(metadataData, 'metadata'));
metadata = metadataData.metadata;
assertEqual(metadata.calculationName, 'PROB');
assertEqual(metadata.quantity, 'physicalDose');
assertEqual(metadata.refScen, 1);
assertEqual(metadata.scenarioDijIx, [1; 2]);
assertEqual(metadata.scenarioCtScenIds, [1; 1]);
assertElementsAlmostEqual(metadata.scenarioWeights, [0.25; 0.75], ...
                          'absolute', 1e-12);
assertTrue(~isempty(dir(fullfile(prob_1.cacheDir, 'scenario_*_voi_*_block_*.mat'))));

function test_probScenarioBatchDiskCacheCleansHashFolderByDefault
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cacheRoot = tempname();
cleanup = onCleanup(@() helper_deleteDirIfExists(cacheRoot));
cfg.CacheRoot = cacheRoot;
cfg.KeepCache = false;

[plnResult, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

dijProb = plnResult.propOpt.dij_prob;
assertFalse(isfield(dijProb, 'cacheDir'));
helper_assertScenarioBatchSizeDisk(dijProb);
assertEqual(numel(helper_listCacheRunDirs(cacheRoot)), 0);

function test_probScenarioBatchMultictMappingUsesCtScenarioIds
[ct, cst, pln, dij, cfg] = helper_multiCtFixture();
pln.propOpt.scen4D = 'all';
[plnResult, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

cfg.SecondPassStrategy = 'recompute';
[plnRecompute, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);

helper_assertProbAlmostEqual(plnResult.propOpt.dij_prob, ...
                             plnRecompute.propOpt.dij_prob);
assertEqual(plnResult.propOpt.dij_prob.scenarioCtScenIds, [1; 2]);

function test_probScenarioBatchSupportsSparseCtScenarioIds
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
ct.numOfCtScen = 3;
ct.refScen = 2;
cfg.refScen = 2;
cst{1, 4} = {[], [1; 2]};
cst{2, 4} = {[], [3; 4]};
scenarioValues = zeros(1, 5);
scenMask = false(3, 1, 1);
scenMask(2, 1, 1) = true;
pln.multScen = helper_fixtureScenarioModel(ct, 2, scenarioValues, ...
                                           [2 1 1], scenMask, 1);
pln.propOpt.scen4D = 2;
dij.physicalDose = cell(3, 1, 1);
dij.physicalDose{2} = sparse([2 0; 0 3; 1 1; 0 4]);

[plnResult, ~] = matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);
dijProb = plnResult.propOpt.dij_prob;

assertEqual(dijProb.refScen, 2);
assertEqual(dijProb.scenarioCtScenIds, 2);
assertElementsAlmostEqual(full(dijProb.expected), ...
                          full(dij.physicalDose{2}), 'absolute', 1e-12);
assertElementsAlmostEqual(full(dijProb.Omega{1}), zeros(2), ...
                          'absolute', 1e-12);

function test_probScenarioBatchRejectsInvalidSecondPassStrategy
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.SecondPassStrategy = 'memory';

assertExceptionThrown(@() matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg), ...
                      'matRad:Error');

function test_probScenarioBatchRejectsInvalidCacheRoot
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.CacheRoot = 1;

assertExceptionThrown(@() matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg), ...
                      'matRad:Error');

function test_probScenarioBatchRejectsCacheRootFile
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.CacheRoot = tempname();
fid = fopen(cfg.CacheRoot, 'w');
fwrite(fid, 'not a folder');
fclose(fid);
cleanup = onCleanup(@() helper_deleteFileIfExists(cfg.CacheRoot));

assertExceptionThrown(@() matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg), ...
                      'matRad:Error');

function test_probRejectsOmegaMemoryAboveLimit
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.MemoryLimitMB = 1e-6;

try
    matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg);
catch ME
    assertEqual(ME.identifier, 'matRad:Error');
    assertTrue(~isempty(strfind(ME.message, 'Probabilistic omega estimated memory')));
    return
end

assertTrue(false);

function test_probRejectsUnresolvedQuantityWithoutLinearField
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.Quantity = 'effect';

assertExceptionThrown(@() matRad_calculateProbabilisticQuantities(ct, cst, [], pln, dij, cfg), ...
                      'matRad:Error');

function [ct, cst, pln, dij, cfg] = helper_singleCtFixture()
ct.numOfCtScen = 1;
ct.cubeDim = [2 2 1];
ct.resolution = struct('x', 1, 'y', 1, 'z', 1);
ct.refScen = 1;

cst = cell(2, 6);
cst = helper_addStructure(cst, 1, 'PTV', 'TARGET', [1; 2]);
cst = helper_addStructure(cst, 2, 'OAR', 'OAR', [3; 4]);

pln.propOpt.quantityOpt = 'physicalDose';
scenarioValues = [0 0 0 0 0; 1 0 0 0 0];
pln.multScen = helper_fixtureScenarioModel(ct, [1; 1], scenarioValues, ...
                                           [1 1 1; 1 2 1], true(1, 2, 1), [0.25; 0.75]);

dij = helper_baseDij(ct.cubeDim, 2);
dij.physicalDose = cell(1, 2, 1);
dij.physicalDose{1} = sparse([1 0; 0 2; 1 1; 0 1]);
dij.physicalDose{2} = sparse([3 0; 0 4; 2 1; 0 5]);
cfg.BatchSize = 2;

function [ct, cst, pln, dij, cfg] = helper_partialSelectionFixture()
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cst{2, 4}{1} = 3;
cst(3, 1:6) = cell(1, 6);
cst = helper_addStructure(cst, 3, 'Spare OAR', 'OAR', 4);
cfg.OARStructSel = 2;

function [ct, cst, pln, dij, cfg] = helper_multiCtFixture()
dim = [2 3 1];
[xGrid, ~, ~] = meshgrid(1:dim(2), 1:dim(1), 1:dim(3));
sourceRows = xGrid(:);

ct.numOfCtScen = 2;
ct.cubeDim = dim;
ct.resolution = struct('x', 1, 'y', 1, 'z', 1);
ct.x = 1:dim(2);
ct.y = 1:dim(1);
ct.z = 1:dim(3);
ct.refScen = 1;
ct.dvfMetadata.dvfType = 'pull';
ct.dvfMetadata.dvfUnits = 'voxel';
ct.dvfMetadata.refScen = 1;
ct.dvfMetadata.referenceCtScen = 1;
ct.dvf = cell(2, 1);
ct.dvf{1} = zeros([3 dim]);
ct.dvf{2} = zeros([3 dim]);
ct.dvf{2}(1, :, :, :) = 1;

cst = cell(2, 6);
cst = helper_addStructure(cst, 1, 'PTV', 'TARGET', [1; 2; 3]);
cst = helper_addStructure(cst, 2, 'OAR', 'OAR', [4; 5; 6]);
cst{1, 4} = cell(1, 2);
cst{2, 4} = cell(1, 2);
cst{1, 4}{1} = [1; 2; 3];
cst{2, 4}{1} = [4; 5; 6];

pln.propOpt.quantityOpt = 'physicalDose';
pln.multScen = helper_fixtureScenarioModel(ct, [1; 2], zeros(2, 5), ...
                                           [1 1 1; 2 1 1], true(2, 1, 1), [0.5; 0.5]);

dij = helper_baseDij(dim, 1);
dij.physicalDose = cell(2, 1, 1);
dij.physicalDose{1} = sparse(sourceRows);
dij.physicalDose{2} = sparse(sourceRows);
cfg.refScen = 1;
cfg.BatchSize = 2;

function [ct, cst, pln, stf] = helper_photonTestDataFixture()
testDataPath = fullfile(fileparts(mfilename('fullpath')), ...
                        '..', 'testData', 'photons_testData.mat');
data = load(testDataPath, 'ct', 'cst', 'pln', 'stf');
ct = data.ct;
cst = data.cst;
pln = data.pln;
stf = data.stf;
pln.propDoseCalc.engine = 'SVDPB';
if ~isfield(pln, 'multScen') || isempty(pln.multScen)
    pln.multScen = matRad_NominalScenario();
end

function scenario = helper_rangeRandomScenario(ct)
scenario = matRad_RandomScenarios(ct);
scenario.nSamples = 2;
scenario.randomSeed = 31;
scenario.scenarioDimensionActive = {'ct', 'range'};

function dij = helper_baseDij(dim, numBixels)
dij.doseGrid.dimensions = dim;
dij.doseGrid.numOfVoxels = prod(dim);
dij.doseGrid.resolution = struct('x', 1, 'y', 1, 'z', 1);
dij.doseGrid.x = 1:dim(2);
dij.doseGrid.y = 1:dim(1);
dij.doseGrid.z = 1:dim(3);
dij.ctGrid = dij.doseGrid;
dij.totalNumOfBixels = numBixels;

function cst = helper_addStructure(cst, rowIx, name, type, voxels)
cst{rowIx, 1} = rowIx;
cst{rowIx, 2} = name;
cst{rowIx, 3} = type;
cst{rowIx, 4} = {voxels(:)};
cst{rowIx, 5} = struct();
cst{rowIx, 6} = {};

function scenarioModel = helper_fixtureScenarioModel(ct, ctScenIds, scenarioValues, ...
                                                     linearMask, scenMask, scenarioWeights)
scenarioModel = matRad_NominalScenario(ct);
dimensions = matRad_createScenarioComponents([1 1 1], 1, 1);
scenForProb = [ctScenIds(:) scenarioValues];
storageSize = size(scenMask);
storageSize(end + 1:size(linearMask, 2)) = 1;
scenarioModel.setScenarioRealizations(dimensions, scenarioValues, ctScenIds, ...
                                      scenarioWeights(:), scenarioWeights(:), scenForProb, linearMask, ...
                                      storageSize, 'legacy-grid');

function omega = helper_manualCenteredOmega(scenarioRows, scenarioWeights, expectedRows, numBixels)
omega = sparse(numBixels, numBixels);
for s = 1:numel(scenarioRows)
    centeredRows = scenarioRows{s} - expectedRows;
    omega = omega + scenarioWeights(s) .* (centeredRows' * centeredRows);
end
omega = sparse(0.5 .* (omega + omega'));

function helper_assertProbAlmostEqual(expectedProb, actualProb)
assertElementsAlmostEqual(full(actualProb.expected), ...
                          full(expectedProb.expected), 'absolute', 1e-12);
assertEqual(size(actualProb.Omega), size(expectedProb.Omega));
for i = 1:numel(expectedProb.Omega)
    if isempty(expectedProb.Omega{i})
        assertTrue(isempty(actualProb.Omega{i}));
    else
        assertElementsAlmostEqual(full(actualProb.Omega{i}), ...
                                  full(expectedProb.Omega{i}), 'absolute', 1e-12);
    end
end
assertEqual(actualProb.voiSubIx, expectedProb.voiSubIx);
assertEqual(actualProb.quantity, expectedProb.quantity);
assertEqual(actualProb.quantityField, expectedProb.quantityField);
assertElementsAlmostEqual(actualProb.scenarioWeights, ...
                          expectedProb.scenarioWeights, 'absolute', 1e-12);

function helper_assertScenarioBatchSizeTotal(dij)
assertTrue(isfield(dij, 'precomputeSize'));
sizeData = dij.precomputeSize;
assertTrue(sizeData.compactBytes > 0);
assertTrue(sizeData.auxiliaryPeakBytes >= 0);
assertElementsAlmostEqual(sizeData.totalPrecomputingBytes, ...
                          sizeData.compactBytes + sizeData.auxiliaryPeakBytes, ...
                          'absolute', 1e-12);

function helper_assertScenarioBatchSizeDisk(dij)
helper_assertScenarioBatchSizeTotal(dij);
sizeData = dij.precomputeSize;
assertTrue(sizeData.diskCachePeakBytes > 0);
assertEqual(sizeData.auxiliaryPeakBytes, ...
            sizeData.diskCachePeakBytes);
assertEqual(sizeData.auxiliaryPeakKind, 'diskCache');
assertEqual(sizeData.secondPassStrategy, 'disk');

function helper_assertScenarioBatchSizeRecompute(dij)
helper_assertScenarioBatchSizeTotal(dij);
sizeData = dij.precomputeSize;
assertTrue(sizeData.memoryTemporaryPeakBytes > 0);
assertEqual(sizeData.auxiliaryPeakBytes, ...
            sizeData.memoryTemporaryPeakBytes);
assertEqual(sizeData.auxiliaryPeakKind, 'memoryTemporary');
assertEqual(sizeData.secondPassStrategy, 'recompute');

function dirs = helper_listCacheRunDirs(cacheRoot)
dirs = {};
if exist(cacheRoot, 'dir') ~= 7
    return
end
listing = dir(cacheRoot);
isRunDir = [listing.isdir] & ~ismember({listing.name}, {'.', '..'});
dirs = {listing(isRunDir).name};

function helper_deleteDirIfExists(path)
if exist(path, 'dir') == 7
    rmdir(path, 's');
end

function helper_deleteFileIfExists(path)
if exist(path, 'file') == 2
    delete(path);
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

function [cleanup, startLogCount] = helper_captureMatRadLog()
matRadCfg = MatRad_Config.instance();
keepLog = matRadCfg.keepLog;
startLogCount = size(matRadCfg.messageLog, 1);
matRadCfg.keepLog = true;
cleanup = onCleanup(@() helper_restoreMatRadLog(matRadCfg, keepLog));

function helper_restoreMatRadLog(matRadCfg, keepLog)
matRadCfg.keepLog = keepLog;

function messages = helper_matRadLogMessages(startLogCount)
if nargin < 1
    startLogCount = 0;
end
matRadCfg = MatRad_Config.instance();
if size(matRadCfg.messageLog, 1) <= startLogCount
    messages = {};
else
    messages = matRadCfg.messageLog(startLogCount + 1:end, 2);
end

function tf = helper_anyMessageContains(messages, needle)
tf = false;
for i = 1:numel(messages)
    if ~isempty(strfind(messages{i}, needle))
        tf = true;
        return
    end
end
