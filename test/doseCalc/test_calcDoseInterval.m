function test_suite = test_calcDoseInterval

test_functions = localfunctions();

initTestSuite;

function test_interval2ScenarioBatchDiskMatchesRecomputeStd
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL2';
[plnDisk, dijDiskContext] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

cfg.SecondPassStrategy = 'recompute';
[plnRecompute, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

diskInterval = plnDisk.propOpt.dij_interval;
recomputeInterval = plnRecompute.propOpt.dij_interval;
assertElementsAlmostEqual(full(diskInterval.center), ...
                          full(recomputeInterval.center), 'absolute', 1e-12);
assertElementsAlmostEqual(full(diskInterval.radius), ...
                          full(recomputeInterval.radius), 'absolute', 1e-12);
assertEqual(diskInterval.precomputeMode, 'scenario-batch');
assertEqual(diskInterval.secondPassStrategy, 'disk');
assertEqual(recomputeInterval.secondPassStrategy, 'recompute');
assertEqual(dijDiskContext.numOfScenarios, 1);

function test_interval2ScenarioBatchPartialSelectionDiskMatchesRecompute
[ct, cst, pln, dij, cfg] = helper_partialSelectionFixture();
cfg.IntervalMode = 'INTERVAL2';
cfg.BatchSize = 1;

[plnDisk, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

assertEqual(plnDisk.propOpt.dij_interval.secondPassStrategy, 'disk');

cfg.SecondPassStrategy = 'recompute';
[plnRecompute, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

assertElementsAlmostEqual(full(plnRecompute.propOpt.dij_interval.center), ...
                          full(plnDisk.propOpt.dij_interval.center), 'absolute', 1e-12);
assertElementsAlmostEqual(full(plnRecompute.propOpt.dij_interval.radius), ...
                          full(plnDisk.propOpt.dij_interval.radius), 'absolute', 1e-12);
assertEqual(plnRecompute.propOpt.dij_interval.secondPassStrategy, ...
            'recompute');

function test_interval2ScenarioBatchSerialDetailedProgressIsAggregated
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL2';
cfg.BatchSize = 1;
cfg.RadiusMode = 'extreme';
cfg.SecondPassStrategy = 'recompute';
cfg.ProgressLevel = 'detailed';

[logCleanup, startLogCount] = helper_captureMatRadLog();
matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

messages = helper_matRadLogMessages(startLogCount);
assertTrue(helper_anyMessageContains(messages, 'Target center progress 4/4'));
assertTrue(helper_anyMessageContains(messages, 'OAR center progress 4/4'));
assertTrue(helper_anyMessageContains(messages, 'Target extreme radius progress 4/4'));
assertTrue(helper_anyMessageContains(messages, 'scenario 2/2, batch 2/2'));

function test_interval2ScenarioBatchUseParallelPrecomputedMatchesSerialStd
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL2';
cfg.BatchSize = 1;

cfg.UseParallel = false;
[plnSerial, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

cfg.UseParallel = true;
[plnParallel, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

assertElementsAlmostEqual(full(plnParallel.propOpt.dij_interval.center), ...
                          full(plnSerial.propOpt.dij_interval.center), 'absolute', 1e-12);
assertElementsAlmostEqual(full(plnParallel.propOpt.dij_interval.radius), ...
                          full(plnSerial.propOpt.dij_interval.radius), 'absolute', 1e-12);
helper_assertScenarioBatchSizeDiskWithoutAuxiliary(plnParallel.propOpt.dij_interval);

function test_interval2ScenarioBatchUseParallelPrecomputedExtremeMatchesSerial
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL2';
cfg.BatchSize = 1;
cfg.RadiusMode = 'extreme';
cfg.SecondPassStrategy = 'recompute';

cfg.UseParallel = false;
[plnSerial, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

cfg.UseParallel = true;
[plnParallel, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

assertElementsAlmostEqual(full(plnParallel.propOpt.dij_interval.center), ...
                          full(plnSerial.propOpt.dij_interval.center), 'absolute', 1e-12);
assertElementsAlmostEqual(full(plnParallel.propOpt.dij_interval.radius), ...
                          full(plnSerial.propOpt.dij_interval.radius), 'absolute', 1e-12);
helper_assertScenarioBatchSizeRecompute(plnParallel.propOpt.dij_interval);

function test_interval2ScenarioBatchRecomputesScenarioDijWithoutDijArgument
[ct, cst, pln, stf] = helper_photonTestDataFixture();
cfg.IntervalMode = 'INTERVAL2';
cfg.SecondPassStrategy = 'disk';
cfg.KeepCache = false;
cfg.BatchSize = 10000;
cfg.targetStructSel = 1;

[plnResult, dijContext] = matRad_calcDoseInterval(ct, cst, stf, pln, cfg);
dijInterval = plnResult.propOpt.dij_interval;

assertEqual(size(dijInterval.center, 1), dijContext.doseGrid.numOfVoxels);
assertEqual(size(dijInterval.center, 2), dijContext.totalNumOfBixels);
assertTrue(nnz(dijInterval.center) > 0);
assertTrue(nnz(dijInterval.radius) > 0);
assertFalse(isfield(dijInterval, 'cacheDir'));
helper_assertScenarioBatchSizeTotal(dijInterval);
assertEqual(dijInterval.intervalMode, 'INTERVAL2');
assertEqual(dijContext.numOfScenarios, 1);

function test_interval2ScenarioBatchCollectTimingReportsRealParallelWhenAvailable
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is unavailable.');
end

[ct, cst, pln, stf] = helper_photonTestDataFixture();
pln.multScen = helper_rangeRandomScenario(ct);
engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
engine.UseParallel = true;
engine.randomSeed = 23;
originalRandomSeed = engine.randomSeed;
originalScenarioIds = engine.multScen.scenarioIds();
pln.propDoseCalc = engine;

cfg.SecondPassStrategy = 'recompute';
cfg.IntervalMode = 'INTERVAL2';
cfg.RadiusMode = 'extreme';
cfg.KeepCache = false;
cfg.BatchSize = 10000;
cfg.targetStructSel = 1;
cfg.UseParallel = true;
cfg.CollectTiming = true;
cfg.ProgressLevel = 'detailed';
cfg.parallelOptions = struct('workerUpperBound', 2);

[logCleanup, startLogCount] = helper_captureMatRadLog();
[plnResult, ~] = matRad_calcDoseInterval(ct, cst, stf, pln, cfg);
timing = plnResult.propOpt.dij_interval.timing;

assertTrue(engine.UseParallel);
assertEqual(engine.randomSeed, originalRandomSeed);
assertEqual(engine.multScen.scenarioIds(), originalScenarioIds);
assertEqual(timing.calculationMode, 'INTERVAL2');
assertEqual(timing.intervalMode, 'INTERVAL2');
if ~timing.parallelScenario.firstPass || ...
        ~timing.parallelScenario.targetExtreme
    moxunit_throw_test_skipped_exception(['Parallel INTERVAL2 ', ...
                                          'scenario-batch could not be activated in this environment.']);
end
pPool = gcp('nocreate');
assertTrue(~isempty(pPool));
assertTrue(pPool.NumWorkers <= 2);
assertTrue(timing.parallelScenario.firstPass);
assertTrue(timing.parallelScenario.targetExtreme);
assertFalse(timing.parallelScenario.oarRadiusFactors);
messages = helper_matRadLogMessages(startLogCount);
assertTrue(helper_anyMessageContains(messages, 'Target center progress'));
assertTrue(helper_anyMessageContains(messages, 'OAR center progress'));
assertTrue(helper_anyMessageContains(messages, 'Target extreme radius progress'));

function test_intervalScenarioBatchRejectsDuplicatePrecomputedDijInputs
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL2';
cfg.PrecomputedDij = dij;

assertExceptionThrown(@() matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg), ...
                      'matRad:Error');

function test_intervalScenarioBatchRejectsUnknownParallelOption
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL2';
cfg.UseParallel = true;
cfg.parallelOptions = struct('workers', 2);

assertExceptionThrown(@() matRad_calcDoseInterval( ...
                                                  ct, cst, [], pln, dij, cfg), 'matRad:Error');

function test_intervalScenarioBatchRejectsFractionalWorkerUpperBound
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL2';
cfg.UseParallel = true;
cfg.parallelOptions = struct('workerUpperBound', 1.5);

assertExceptionThrown(@() matRad_calcDoseInterval( ...
                                                  ct, cst, [], pln, dij, cfg), 'matRad:Error');

function test_interval2ScenarioBatchExtremeDiskMatchesRecompute
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL2';
cfg.RadiusMode = 'extreme';
[plnDisk, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

cfg.SecondPassStrategy = 'recompute';
[plnResult, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

assertElementsAlmostEqual(full(plnResult.propOpt.dij_interval.center), ...
                          full(plnDisk.propOpt.dij_interval.center), 'absolute', 1e-12);
assertElementsAlmostEqual(full(plnResult.propOpt.dij_interval.radius), ...
                          full(plnDisk.propOpt.dij_interval.radius), 'absolute', 1e-12);
helper_assertScenarioBatchSizeRecompute(plnResult.propOpt.dij_interval);

function test_intervalTargetExtremeParallelFallsBackToSerialWhenOnlyParallelEstimateExceedsLimit
ctx.scenarioDijIx = (1:4)';
ctx.targetRows = (1:1000)';
ctx.numBixels = 1000;

batch = struct();
batch.rows = 1;
targetBatches = {batch};

cfg.MemoryLimitMB = 2;

[logCleanup, startLogCount] = helper_captureMatRadLog();
useParallel = DoseInterval.matRad_guardDoseIntervalTargetExtremeMemory(ctx, ...
                                                                       targetBatches, cfg, true, ...
                                                                       MatRad_Config.instance());

messages = helper_matRadLogMessages(startLogCount);
assertFalse(useParallel);
assertTrue(helper_anyMessageContains(messages, ...
                                     'Falling back to serial target extreme radius'));

function test_intervalTargetExtremeRejectsWhenSerialEstimateExceedsLimit
ctx.scenarioDijIx = (1:2)';
ctx.targetRows = (1:1000)';
ctx.numBixels = 1000;

batch = struct();
batch.rows = 1;
targetBatches = {batch};

cfg.MemoryLimitMB = 1e-6;

try
    DoseInterval.matRad_guardDoseIntervalTargetExtremeMemory(ctx, targetBatches, cfg, ...
                                                            false, MatRad_Config.instance());
catch ME
    assertEqual(ME.identifier, 'matRad:Error');
    assertTrue(~isempty(strfind(ME.message, ...
                                'INTERVAL target extreme radius estimated memory')));
    return
end

assertTrue(false);

function test_interval3ScenarioBatchUseParallelPrecomputedMatchesSerial
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL3';
cfg.KMode = 'static';
cfg.KMax = 2;
cfg.BatchSize = 1;

cfg.UseParallel = false;
[plnSerial, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

cfg.UseParallel = true;
[plnParallel, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

serialInterval = plnSerial.propOpt.dij_interval;
parallelInterval = plnParallel.propOpt.dij_interval;
assertElementsAlmostEqual(full(parallelInterval.center), ...
                          full(serialInterval.center), 'absolute', 1e-12);
assertElementsAlmostEqual(full(parallelInterval.radius), ...
                          full(serialInterval.radius), 'absolute', 1e-12);
assertElementsAlmostEqual(parallelInterval.OARRadiusRank, ...
                          serialInterval.OARRadiusRank, 'absolute', 1e-12);
for i = 1:numel(serialInterval.OARRadiusFactor)
    parallelCovariance = full(helper_reconstructOARRadiusCovariance(parallelInterval, i));
    serialCovariance = full(helper_reconstructOARRadiusCovariance(serialInterval, i));
    assertElementsAlmostEqual(parallelCovariance, serialCovariance, 'absolute', 1e-12);
end

function test_interval3ScenarioBatchDiskCacheKeepsDistinctHashFolders
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.RadiusMode = 'extreme';
cfg.KMode = 'static';
cfg.KMax = 2;

cacheRoot = tempname();
cleanup = onCleanup(@() helper_deleteDirIfExists(cacheRoot));
cfg.IntervalMode = 'INTERVAL3';
cfg.CacheRoot = cacheRoot;
cfg.KeepCache = true;

[plnRun1, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);
[plnRun2, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

interval1 = plnRun1.propOpt.dij_interval;
interval2 = plnRun2.propOpt.dij_interval;
assertTrue(isfield(interval1, 'cacheDir'));
assertTrue(isfield(interval2, 'cacheDir'));
helper_assertScenarioBatchSizeDisk(interval1);
helper_assertScenarioBatchSizeDisk(interval2);
assertFalse(strcmp(interval1.cacheDir, interval2.cacheDir));
assertEqual(exist(interval1.cacheDir, 'dir'), 7);
assertEqual(exist(interval2.cacheDir, 'dir'), 7);
assertEqual(exist(fullfile(interval1.cacheDir, 'metadata.mat'), 'file'), 2);
metadataData = load(fullfile(interval1.cacheDir, 'metadata.mat'));
assertTrue(isfield(metadataData, 'metadata'));
metadata = metadataData.metadata;
assertEqual(metadata.calculationName, 'dose interval');
assertEqual(metadata.intervalMode, 'INTERVAL3');
assertEqual(metadata.radiusMode, 'extreme');
assertEqual(metadata.quantity, 'physicalDose');
assertEqual(metadata.refScen, 1);
assertEqual(metadata.scenarioDijIx, [1; 2]);
assertEqual(metadata.scenarioCtScenIds, [1; 1]);
assertElementsAlmostEqual(metadata.scenarioWeights, [0.25; 0.75], ...
                          'absolute', 1e-12);
assertTrue(~isempty(dir(fullfile(interval1.cacheDir, 'scenario_*_oar_block_*.mat'))));

assertElementsAlmostEqual(full(interval1.center), ...
                          full(interval2.center), 'absolute', 1e-12);
assertElementsAlmostEqual(full(interval1.radius), ...
                          full(interval2.radius), 'absolute', 1e-12);
assertElementsAlmostEqual(full(helper_reconstructOARRadiusCovariance(interval1, 1)), ...
                          full(helper_reconstructOARRadiusCovariance(interval2, 1)), 'absolute', 1e-12);

function test_interval3ScenarioBatchDiskCacheCleansHashFolderByDefault
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cacheRoot = tempname();
cleanup = onCleanup(@() helper_deleteDirIfExists(cacheRoot));
cfg.IntervalMode = 'INTERVAL3';
cfg.CacheRoot = cacheRoot;
cfg.KeepCache = false;
cfg.KMode = 'static';
cfg.KMax = 2;

[plnResult, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

dijInterval = plnResult.propOpt.dij_interval;
assertFalse(isfield(dijInterval, 'cacheDir'));
helper_assertScenarioBatchSizeDisk(dijInterval);
assertEqual(numel(helper_listCacheRunDirs(cacheRoot)), 0);

function test_interval2ScenarioBatchMultictMappingUsesCtScenarioIds
[ct, cst, pln, dij, cfg, expectedCenter] = helper_multiCtFixture(1);
pln.propOpt.scen4D = 'all';
cfg.IntervalMode = 'INTERVAL2';

[plnResult, ~] = matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg);

assertElementsAlmostEqual(full(plnResult.propOpt.dij_interval.center(:, 1)), ...
                          expectedCenter, 'absolute', 1e-12);
assertEqual(plnResult.propOpt.dij_interval.scenarioCtScenIds, [1; 2]);

function test_intervalScenarioBatchRequiresIntervalMode
[ct, cst, pln, ~, cfg] = helper_singleCtFixture();

assertExceptionThrown(@() matRad_calcDoseInterval(ct, cst, [], pln, cfg), ...
                      'matRad:Error');

function test_intervalScenarioBatchRejectsInvalidSecondPassStrategy
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL2';
cfg.SecondPassStrategy = 'memory';

assertExceptionThrown(@() matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg), ...
                      'matRad:Error');

function test_intervalScenarioBatchRejectsInvalidCacheRoot
[ct, cst, pln, ~, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL2';
cfg.CacheRoot = 1;

assertExceptionThrown(@() matRad_calcDoseInterval(ct, cst, [], pln, cfg), ...
                      'matRad:Error');

function test_intervalScenarioBatchRejectsCacheRootFile
[ct, cst, pln, dij, cfg] = helper_singleCtFixture();
cfg.IntervalMode = 'INTERVAL3';
cfg.CacheRoot = tempname();
fid = fopen(cfg.CacheRoot, 'w');
fwrite(fid, 'not a folder');
fclose(fid);
cleanup = onCleanup(@() helper_deleteFileIfExists(cfg.CacheRoot));

assertExceptionThrown(@() matRad_calcDoseInterval(ct, cst, [], pln, dij, cfg), ...
                      'matRad:Error');

function covariance = helper_reconstructOARRadiusCovariance(dijInterval, intervalIx)
factor = dijInterval.OARRadiusFactor{intervalIx};
covariance = factor * factor';

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

function [ct, cst, pln, dij, cfg, expectedCenter] = helper_multiCtFixture(refScen)
dim = [2 3 1];
[xGrid, ~, ~] = meshgrid(1:dim(2), 1:dim(1), 1:dim(3));
sourceRows = xGrid(:);
mappedRows = zeros(dim);
mappedRows(:, 2:end, :) = xGrid(:, 1:end - 1, :);
mappedRows = mappedRows(:);

ct.numOfCtScen = refScen + 1;
ct.cubeDim = dim;
ct.resolution = struct('x', 1, 'y', 1, 'z', 1);
ct.x = 1:dim(2);
ct.y = 1:dim(1);
ct.z = 1:dim(3);
ct.refScen = refScen;
ct.dvfMetadata.dvfType = 'pull';
ct.dvfMetadata.dvfUnits = 'voxel';
ct.dvfMetadata.refScen = refScen;
ct.dvfMetadata.referenceCtScen = refScen;
ct.dvf = cell(ct.numOfCtScen, 1);
ct.dvf{refScen} = zeros([3 dim]);
ct.dvf{refScen + 1} = zeros([3 dim]);
ct.dvf{refScen + 1}(1, :, :, :) = 1;

cst = cell(1, 6);
cst = helper_addStructure(cst, 1, 'PTV', 'TARGET', (1:prod(dim))');
cst{1, 4} = cell(1, ct.numOfCtScen);
cst{1, 4}{refScen} = (1:prod(dim))';

pln.propOpt.quantityOpt = 'physicalDose';
scenMask = false(ct.numOfCtScen, 1, 1);
scenMask(refScen) = true;
scenMask(refScen + 1) = true;
pln.multScen = helper_fixtureScenarioModel(ct, [refScen; refScen + 1], zeros(2, 5), ...
                                           [refScen 1 1; refScen + 1 1 1], scenMask, [0.5; 0.5]);

dij = helper_baseDij(dim, 1);
dij.physicalDose = cell(ct.numOfCtScen, 1, 1);
dij.physicalDose{refScen} = sparse(sourceRows);
dij.physicalDose{refScen + 1} = sparse(sourceRows);
cfg.refScen = refScen;
cfg.BatchSize = 2;

expectedCenter = 0.5 * sourceRows + 0.5 * mappedRows;

function [ct, cst, pln, dij, cfg, expectedCenter] = helper_multiCtYDisplacementFixture()
dim = [3 2 1];
[~, yGrid, ~] = meshgrid(1:dim(2), 1:dim(1), 1:dim(3));
sourceRows = yGrid(:);
mappedRows = zeros(dim);
mappedRows(2:end, :, :) = yGrid(1:end - 1, :, :);
mappedRows = mappedRows(:);

[ct, cst, pln, dij, cfg] = helper_multiCtMappingBaseFixture(dim, sourceRows);
ct.dvf{2}(2, :, :, :) = 1;

expectedCenter = 0.5 * sourceRows + 0.5 * mappedRows;

function [ct, cst, pln, dij, cfg, expectedCenter] = helper_multiCtMillimeterDvfFixture()
dim = [2 3 1];
[xGrid, ~, ~] = meshgrid(1:dim(2), 1:dim(1), 1:dim(3));
sourceRows = xGrid(:);
mappedRows = zeros(dim);
mappedRows(:, 2:end, :) = xGrid(:, 1:end - 1, :);
mappedRows = mappedRows(:);

[ct, cst, pln, dij, cfg] = helper_multiCtMappingBaseFixture(dim, sourceRows);
ct.resolution = struct('x', 2, 'y', 1, 'z', 1);
ct.dvfMetadata.dvfUnits = 'mm';
ct.dvf{2}(1, :, :, :) = 2;
dij.doseGrid.resolution = struct('x', 2, 'y', 1, 'z', 1);
dij.ctGrid = dij.doseGrid;

expectedCenter = 0.5 * sourceRows + 0.5 * mappedRows;

function [ct, cst, pln, dij, cfg] = helper_multiCtMappingBaseFixture(dim, sourceRows)
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

cst = cell(1, 6);
cst = helper_addStructure(cst, 1, 'PTV', 'TARGET', (1:prod(dim))');
cst{1, 4} = cell(1, 2);
cst{1, 4}{1} = (1:prod(dim))';

pln.propOpt.quantityOpt = 'physicalDose';
pln.multScen = helper_fixtureScenarioModel(ct, [1; 2], zeros(2, 5), ...
                                           [1 1 1; 2 1 1], true(2, 1, 1), [0.5; 0.5]);

dij = helper_baseDij(dim, 1);
dij.physicalDose = cell(2, 1, 1);
dij.physicalDose{1} = sparse(sourceRows);
dij.physicalDose{2} = sparse(sourceRows);
cfg.refScen = 1;
cfg.BatchSize = 2;

function [ct, cst, pln, dij, cfg, expectedCenter] = helper_multiCtResamplingFixture()
ctDim = [2 4 2];
doseDim = [2 2 2];

ct.numOfCtScen = 2;
ct.cubeDim = ctDim;
ct.resolution = struct('x', 1, 'y', 1, 'z', 1);
ct.refScen = 1;
ct.dvfMetadata.dvfType = 'pull';
ct.dvfMetadata.dvfUnits = 'mm';
ct.dvfMetadata.refScen = 1;
ct.dvfMetadata.referenceCtScen = 1;
ct.dvf = cell(2, 1);
ct.dvf{1} = zeros([3 ctDim]);
ct.dvf{2} = zeros([3 ctDim]);

cst = cell(1, 6);
cst = helper_addStructure(cst, 1, 'PTV', 'TARGET', (1:prod(ctDim))');

pln.propOpt.quantityOpt = 'physicalDose';
pln.multScen = helper_fixtureScenarioModel(ct, [1; 2], zeros(2, 5), ...
                                           [1 1 1; 2 1 1], true(2, 1, 1), [0.5; 0.5]);

dij = helper_baseDij(doseDim, 1);
dij.ctGrid.dimensions = ctDim;
dij.ctGrid.numOfVoxels = prod(ctDim);
dij.ctGrid.resolution = ct.resolution;
dij.ctGrid.x = 1:ctDim(2);
dij.ctGrid.y = 1:ctDim(1);
dij.ctGrid.z = 1:ctDim(3);
dij.doseGrid.resolution = struct('x', 2, 'y', 1, 'z', 1);
dij.doseGrid.x = 1:2:ctDim(2);
dij.doseGrid.y = 1:ctDim(1);
dij.doseGrid.z = 1:ctDim(3);
dij.physicalDose = cell(2, 1, 1);
dij.physicalDose{1} = sparse((1:prod(doseDim))');
dij.physicalDose{2} = sparse((prod(doseDim) + 1:2 * prod(doseDim))');
cfg.refScen = 1;
cfg.BatchSize = 2;
expectedCenter = 0.5 * ((1:prod(doseDim))' + (prod(doseDim) + 1:2 * prod(doseDim))');

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

function helper_assertScenarioBatchSizeDiskWithoutAuxiliary(dij)
helper_assertScenarioBatchSizeTotal(dij);
sizeData = dij.precomputeSize;
assertEqual(sizeData.diskCachePeakBytes, 0);
assertEqual(sizeData.memoryTemporaryPeakBytes, 0);
assertEqual(sizeData.auxiliaryPeakBytes, 0);
assertElementsAlmostEqual(sizeData.totalPrecomputingBytes, ...
                          sizeData.compactBytes, 'absolute', 1e-12);
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
