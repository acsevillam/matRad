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
