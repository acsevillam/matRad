function test_suite = test_estimateSamplingMemory

test_functions = localfunctions();

initTestSuite;

function test_context_struct_returns_expected_memory_fields
    samplingContext = samplingMemoryFixture();

    samplingMemoryEstimate = matRad_estimateSamplingMemory(samplingContext);

    assertEqual(samplingMemoryEstimate.estimateBasis,'nominalScenarioProxy');
    assertEqual(samplingMemoryEstimate.numSamples,samplingContext.numSamples);
    assertEqual(samplingMemoryEstimate.numVoxels,numel(samplingContext.subIx));
    assertEqual(samplingMemoryEstimate.sampleDoseBytes,numel(samplingContext.subIx) * 4);
    assertEqual(samplingMemoryEstimate.sampleDoseStorageBytes, ...
        samplingContext.samplingDoseStorageBytes);
    assertTrue(samplingMemoryEstimate.rawWorkerBytes > 0);
    assertEqual(samplingMemoryEstimate.workerLimit,[]);

function test_requires_context_struct
    assertExceptionThrown(@() matRad_estimateSamplingMemory([]),'matRad:Error');

function test_requires_expected_context_fields
    samplingContext = samplingMemoryFixture();
    samplingContext = rmfield(samplingContext,'subIx');

    assertExceptionThrown(@() matRad_estimateSamplingMemory(samplingContext),'matRad:Error');

function test_worker_upper_bound_limits_memory_limited_worker_count
    [maxWorkers,memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount(1, ...
        'numTasks',100, ...
        'workerUpperBound',3, ...
        'limitToDefaultPool',false);

    if isempty(maxWorkers)
        moxunit_throw_test_skipped_exception('System memory information is unavailable.');
    end

    assertEqual(maxWorkers,3);
    assertEqual(memoryEstimate.workerUpperBound,3);
    assertEqual(memoryEstimate.maxWorkers,3);

function test_worker_upper_bound_one_limits_to_serial_worker_count
    [maxWorkers,memoryEstimate] = matRad_estimateMemoryLimitedWorkerCount(1, ...
        'numTasks',100, ...
        'workerUpperBound',1, ...
        'limitToDefaultPool',false);

    if isempty(maxWorkers)
        moxunit_throw_test_skipped_exception('System memory information is unavailable.');
    end

    assertEqual(maxWorkers,1);
    assertEqual(memoryEstimate.workerUpperBound,1);
    assertEqual(memoryEstimate.maxWorkers,1);

function test_worker_upper_bound_rejects_non_integer
    assertExceptionThrown(@() matRad_estimateMemoryLimitedWorkerCount(1, ...
        'workerUpperBound',1.5),'matRad:Error');

function test_worker_upper_bound_rejects_zero
    assertExceptionThrown(@() matRad_estimateMemoryLimitedWorkerCount(1, ...
        'workerUpperBound',0),'matRad:Error');

function test_worker_upper_bound_rejects_inf
    assertExceptionThrown(@() matRad_estimateMemoryLimitedWorkerCount(1, ...
        'workerUpperBound',Inf),'matRad:Error');

function test_worker_upper_bound_rejects_text
    assertExceptionThrown(@() matRad_estimateMemoryLimitedWorkerCount(1, ...
        'workerUpperBound','2'),'matRad:Error');

function samplingContext = samplingMemoryFixture()
    samplingContext = struct();
    samplingContext.ct = struct('cubeDim',[2 2 1]);
    samplingContext.stf = struct('ray',1);
    samplingContext.cst = cell(1,6);
    samplingContext.cstEval = cell(1,6);
    samplingContext.pln = struct();
    samplingContext.pln.bioParam = struct('model','none');
    samplingContext.pln.multScen = struct();
    samplingContext.pln.multScen.extractSingleScenario = @extractSingleScenario;
    samplingContext.w = [1; 2];
    samplingContext.subIx = [1; 3];
    samplingContext.samplingCtScenIds = [1; 1];
    samplingContext.dvhPoints = [0 1 2];
    samplingContext.refGy = [0 2];
    samplingContext.refVol = [2 50 98];
    samplingContext.resultGUInomScen = struct('physicalDose',ones([2 2 1]));
    samplingContext.doseMapping = struct('enabled',false,'method','none');
    samplingContext.refScen = 1;
    samplingContext.samplingDoseStorageBytes = 32;
    samplingContext.numSamples = 2;

function scenario = extractSingleScenario(~)
    scenario = struct();
    scenario.getRangeShift = @(~) [0 0];
    scenario.getSetupShift = @(~) [0 0 0];
