function test_suite = test_samplingMemoryCalibration

test_functions = localfunctions();

initTestSuite;

function test_undermeasuredCalibrationKeepsStaticEstimate
staticWorkerBytes = 20 * 1024^3;
measuredWorkerBytes = 116 * 1024^2;
[memoryEstimate, calibration, resourceConfig] = ...
    helper_createCalibrationFixture(staticWorkerBytes, measuredWorkerBytes);

[memoryEstimate, calibration] = ...
    matRad_applySamplingMemoryCalibration(memoryEstimate, calibration, resourceConfig);

assertEqual(memoryEstimate.rawWorkerBytes, staticWorkerBytes);
assertEqual(calibration.action, 'undermeasured');
assertFalse(calibration.usedForPlanning);

function test_highCalibrationRaisesStaticEstimate
staticWorkerBytes = 20 * 1024^3;
measuredWorkerBytes = 24 * 1024^3;
[memoryEstimate, calibration, resourceConfig] = ...
    helper_createCalibrationFixture(staticWorkerBytes, measuredWorkerBytes);

[memoryEstimate, calibration] = ...
    matRad_applySamplingMemoryCalibration(memoryEstimate, calibration, resourceConfig);

assertEqual(memoryEstimate.rawWorkerBytes, measuredWorkerBytes);
assertEqual(calibration.action, 'raised');
assertTrue(calibration.usedForPlanning);

function test_reliableModerateCalibrationLowersToCalibratedFloor
staticWorkerBytes = 8 * 1024^3;
measuredWorkerBytes = 3 * 1024^3;
[memoryEstimate, calibration, resourceConfig] = ...
    helper_createCalibrationFixture(staticWorkerBytes, measuredWorkerBytes);

[memoryEstimate, calibration] = ...
    matRad_applySamplingMemoryCalibration(memoryEstimate, calibration, resourceConfig);

assertEqual(memoryEstimate.rawWorkerBytes, 4 * 1024^3);
assertEqual(calibration.action, 'lowered');
assertTrue(calibration.usedForPlanning);

function test_unreliableCalibrationKeepsStaticEstimate
staticWorkerBytes = 12 * 1024^3;
[memoryEstimate, calibration, resourceConfig] = ...
    helper_createCalibrationFixture(staticWorkerBytes, 20 * 1024^3);
calibration.measurementReliable = false;

[memoryEstimate, calibration] = ...
    matRad_applySamplingMemoryCalibration(memoryEstimate, calibration, resourceConfig);

assertEqual(memoryEstimate.rawWorkerBytes, staticWorkerBytes);
assertEqual(calibration.action, 'unreliable');
assertFalse(calibration.usedForPlanning);

function test_samplingParallelProgressDoesNotUseParforProgress
matRad_cfg = MatRad_Config.instance();
samplingPath = fullfile(matRad_cfg.matRadRoot, 'matRad', 'planAnalysis', ...
                        'matRad_sampling.m');
samplingText = fileread(samplingPath);

assertTrue(isempty(strfind(samplingText, 'parfor_progress')));

function [memoryEstimate, calibration, resourceConfig] = ...
    helper_createCalibrationFixture(staticWorkerBytes, measuredWorkerBytes)

memoryEstimate = struct();
memoryEstimate.rawWorkerBytes = staticWorkerBytes;
memoryEstimate.estimateBasis = 'samplingForwardDoseStatic';

calibration = struct();
calibration.measurementReliable = true;
calibration.measuredWorkerBytes = measuredWorkerBytes;
calibration.staticWorkerBytes = staticWorkerBytes;
calibration.calibratedWorkerBytes = [];
calibration.calibratedMinWorkerBytes = [];
calibration.minReliableMeasuredBytes = [];
calibration.minReliableReductionRatio = [];
calibration.allowReduction = true;
calibration.usedForPlanning = false;
calibration.action = 'kept';
calibration.undermeasurementReason = '';

resourceConfig = struct();
resourceConfig.allowCalibrationToReduceWorkerMemory = true;
resourceConfig.calibratedMinForwardDoseWorkerMemoryBytes = 4 * 1024^3;
resourceConfig.calibrationMinReliableMeasuredBytes = 1 * 1024^3;
resourceConfig.calibrationMinReliableReductionRatio = 0.50;
