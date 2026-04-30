function test_suite = test_calcSamplingDoseStatistics

test_functions = localfunctions();

initTestSuite;

function test_partial_sampling_sets_mask_and_nan_values
    ct = struct();
    ct.cubeDim = [3 2 1];

    pln = struct();
    pln.subIx = [1;4];

    mSampDose = single([1 3;2 4]);
    scenWeights = [1;3];

    doseStat = matRad_calcSamplingDoseStatistics(ct,pln,mSampDose,scenWeights);

    expectedMask = false(ct.cubeDim);
    expectedMask(pln.subIx) = true;

    assertEqual(doseStat.sampleMask,expectedMask);
    assertEqual(doseStat.sampledVoxelIndices,pln.subIx);
    assertElementsAlmostEqual(doseStat.sampleCoverageFraction,2/6);
    assertElementsAlmostEqual(doseStat.meanCube(1),2);
    assertElementsAlmostEqual(doseStat.meanCubeW(1),2.5);
    assertTrue(isnan(doseStat.meanCube(2)));
    assertTrue(isnan(doseStat.stdCubeW(2)));

function test_rejects_mismatched_sample_rows
    ct = struct();
    ct.cubeDim = [3 1 1];

    pln = struct();
    pln.subIx = [1;2];

    assertExceptionThrown(@() matRad_calcSamplingDoseStatistics(ct,pln,single(zeros(1,2)),[1;1]), ...
        'matRad:Error');

function test_rejects_mismatched_scenario_weights
    ct = struct();
    ct.cubeDim = [3 1 1];

    pln = struct();
    pln.subIx = [1;2];

    assertExceptionThrown(@() matRad_calcSamplingDoseStatistics(ct,pln,single(zeros(2,2)),1), ...
        'matRad:Error');
