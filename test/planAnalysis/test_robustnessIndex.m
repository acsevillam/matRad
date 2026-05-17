function test_suite = test_robustnessIndex

test_functions = localfunctions();

initTestSuite;

function test_index1UsesMeanAndStdCriteria
ct = helper_createCt();
cst = helper_createTargetCst([1; 2; 3; 4], 5);
meanCube = reshape([5; 5.25; 5; 6], ct.cubeDim);
stdCube = reshape([0; 0; 0.125; 0.5], ct.cubeDim);

[robustnessCube, robPassRate] = matRad_robustnessIndex(meanCube, stdCube, ...
                                                       5, [10 5], ...
                                                       'index1', ct, cst);

expectedCube = reshape([0; 0.5; 0.5; sqrt(8)], ct.cubeDim);
assertElementsAlmostEqual(robustnessCube, expectedCube, 'absolute', 1e-12);
assertElementsAlmostEqual(robPassRate, 75, 'absolute', 1e-12);

function test_plotRobustnessIndexShowsMetricTextbox
figureCleaner = onCleanup(@() close('all'));
ct = helper_createCt();
cst = helper_createTargetCst([1; 2; 3; 4], 5);
meanCube = reshape([5; 5.25; 5; 6], ct.cubeDim);
stdCube = reshape([0; 0; 0.125; 0.5], ct.cubeDim);

[~, ~, robustnessFig] = matRad_robustnessIndex(meanCube, stdCube, ...
                                               5, [10 5], ...
                                               'index1', ct, cst, 1);

textHandle = findobj(robustnessFig, 'Tag', 'matRadMetricTextbox');
assertEqual(numel(textHandle), 1);
assertEqual(get(textHandle, 'String'), 'RI = 0.7500');

function test_samplingRobustnessPlotWithoutAxesHandleCreatesIndependentFigure
figureCleaner = onCleanup(@() close('all'));
ct = helper_createPlotCt();
cst = cell(0, 6);
analysis.index1.robustnessCube = zeros(ct.cubeDim);
analysis.index1.robustnessIndex = 1;
existingFig = figure('Visible', 'off');
axes('Parent', existingFig);

newFig = matRad_plotSamplingRobustnessAnalysis(analysis, ct, cst, 1, ...
                                               'method', 'index1');

assertFalse(isequal(newFig, existingFig));
plotAxes = findobj(newFig, 'Type', 'Axes');
assertFalse(isempty(plotAxes));

function test_index2UsesBinaryCriteria
ct = helper_createCt();
cst = helper_createTargetCst([1; 2; 3; 4], 5);
meanCube = reshape([5; 5.25; 5; 6], ct.cubeDim);
stdCube = reshape([0; 0; 0.125; 0.5], ct.cubeDim);

[robustnessCube, robPassRate] = matRad_robustnessIndex(meanCube, stdCube, ...
                                                       5, [10 5], ...
                                                       'index2', ct, cst);

expectedCube = reshape([true; true; true; false], ct.cubeDim);
assertEqual(robustnessCube, expectedCube);
assertElementsAlmostEqual(robPassRate, 75, 'absolute', 1e-12);

function test_samplingRobustnessSkipsPartialTargets
ct = helper_createCt();
cst = helper_createTargetCst([1; 4], 10);
pln = helper_createPlan(1);
meanCube = 10 * ones(ct.cubeDim);
stdCube = zeros(ct.cubeDim);
sampleMask = true(ct.cubeDim);
sampleMask(4) = false;

analysis = matRad_samplingRobustnessAnalysis(meanCube, stdCube, [5 5], ...
                                             ct, cst, pln, [], ...
                                             'sampleMask', sampleMask);

assertFalse(analysis.targets(1).isEvaluable);
assertEqual(analysis.targets(1).samplingStatus, 'partialSamplingCoverage');
assertEqual(analysis.targets(1).numUnsampledVoxels, 1);
assertTrue(isempty(analysis.index1.robPassRate));
assertTrue(all(isnan(analysis.index1.robustnessCube(:))));

function test_samplingRobustnessUsesPerFractionTargetReferenceDose
ct = helper_createCt();
cst = helper_createTargetCst([1; 2; 3; 4], 10);
pln = helper_createPlan(2);
meanCube = 5 * ones(ct.cubeDim);
stdCube = zeros(ct.cubeDim);

analysis = matRad_samplingRobustnessAnalysis(meanCube, stdCube, [1 1], ...
                                             ct, cst, pln);

assertElementsAlmostEqual(analysis.targets(1).refDose, 5, 'absolute', 1e-12);
assertTrue(analysis.targets(1).isEvaluable);
assertElementsAlmostEqual(analysis.index1.robPassRate, 100, 'absolute', 1e-12);
assertElementsAlmostEqual(analysis.index2.robPassRate, 100, 'absolute', 1e-12);

function ct = helper_createCt()
ct.cubeDim = [2 2];
ct.refScen = 1;
ct.resolution.x = 1;
ct.resolution.y = 1;
ct.resolution.z = 1;

function ct = helper_createPlotCt()
ct.cubeDim = [2 2 2];
ct.refScen = 1;
ct.resolution.x = 1;
ct.resolution.y = 1;
ct.resolution.z = 1;
ct.cubeHU = {zeros(ct.cubeDim)};

function cst = helper_createTargetCst(voxels, refDoseTotal)
cst = cell(1, 6);
cst{1, 1} = 1;
cst{1, 2} = 'target';
cst{1, 3} = 'TARGET';
cst{1, 4} = {voxels(:)};
cst{1, 5} = struct('Visible', true);
cst{1, 6} = {DoseObjectives.matRad_SquaredDeviation(1, refDoseTotal)};

function pln = helper_createPlan(numOfFractions)
pln.numOfFractions = numOfFractions;
