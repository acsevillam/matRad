function test_suite = test_expectedDoseDifferenceAnalysis

test_functions = localfunctions();

initTestSuite;

function test_weightedOverUnderAndNearReferenceProbabilities
referenceCube = reshape([1; 1; 2; 4], [2 2]);
sampleDoseMatrix = [2.00 0.00 2.00; ...
                    1.00 2.00 0.00; ...
                    1.95 2.00 2.05; ...
                    4.00 4.20 3.70];
weights = [0.2; 0.3; 0.5];

analysis = matRad_expectedDoseDifferenceAnalysis(sampleDoseMatrix, referenceCube, ...
                                                 'scenWeights', weights, ...
                                                 'tolerance', 0.1, ...
                                                 'quantity', 'physicalDose', ...
                                                 'evaluationModeBase', 'perFraction', ...
                                                 'referenceName', 'nominalDose');

assertEqual(analysis.status, 'computedFullCube');
assertEqual(analysis.quantity, 'physicalDose');
assertEqual(analysis.evaluationModeBase, 'perFraction');
assertEqual(analysis.referenceName, 'nominalDose');
assertElementsAlmostEqual(analysis.overReferenceProbabilityCube(1), 0.7, 'absolute', 1e-12);
assertElementsAlmostEqual(analysis.underReferenceProbabilityCube(1), 0.3, 'absolute', 1e-12);
assertElementsAlmostEqual(analysis.nearReferenceProbabilityCube(3), 1.0, 'absolute', 1e-12);
assertElementsAlmostEqual(analysis.signedReferenceProbabilityCube(4), -0.2, 'absolute', 1e-12);
assertElementsAlmostEqual(analysis.overReferenceExpectedDoseDifferenceCube(1), 0.7, ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(analysis.underReferenceExpectedDoseDifferenceCube(1), 0.3, ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(analysis.signedExpectedDoseDifferenceCube(1), 0.4, ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(analysis.signedExpectedDoseDifferenceCube(4), -0.09, ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(analysis.summary.maxOverReferenceProbability, 0.7, 'absolute', 1e-12);
assertElementsAlmostEqual(analysis.summary.maxUnderReferenceProbability, 0.5, 'absolute', 1e-12);
assertElementsAlmostEqual(analysis.summary.maxOverReferenceExpectedDoseDifference, 0.7, ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(analysis.summary.maxAbsSignedExpectedDoseDifference, 0.4, ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(analysis.doseWindow, [-0.4 0.4], 'absolute', 1e-12);

function test_configuredDoseWindowIsStoredInAnalysis
referenceCube = ones(2, 2);
sampleDoseMatrix = [1.4 0.9; ...
                    1.1 1.2; ...
                    0.8 1.0; ...
                    1.0 1.0];

analysis = matRad_expectedDoseDifferenceAnalysis(sampleDoseMatrix, referenceCube, ...
                                                 'doseWindow', [-2 3]);

assertElementsAlmostEqual(analysis.doseWindow, [-2 3], 'absolute', 1e-12);

function test_partialMaskLeavesUnsampledVoxelsAsNaN
referenceCube = reshape([1; 1; 2; 4], [2 2]);
sampleDoseMatrix = [2 0 2; ...
                    4 5 3];
weights = [0.2; 0.3; 0.5];
sampleMask = false(size(referenceCube));
sampleMask([1 4]) = true;

analysis = matRad_expectedDoseDifferenceAnalysis(sampleDoseMatrix, referenceCube, ...
                                                 'sampledVoxelIndices', [1; 4], ...
                                                 'sampleMask', sampleMask, ...
                                                 'scenWeights', weights);

assertEqual(analysis.status, 'computedPartialMask');
assertElementsAlmostEqual(analysis.sampleCoverageFraction, 0.5, 'absolute', 1e-12);
assertElementsAlmostEqual(analysis.overReferenceProbabilityCube(1), 0.7, 'absolute', 1e-12);
assertElementsAlmostEqual(analysis.underReferenceProbabilityCube(4), 0.5, 'absolute', 1e-12);
assertTrue(isnan(analysis.overReferenceProbabilityCube(2)));
assertTrue(isnan(analysis.underReferenceProbabilityCube(3)));

function test_plotUsesSignedExpectedDoseDifference
figureCleaner = onCleanup(@() close('all'));
referenceCube = reshape([1; 1; 2; 4], [2 2]);
sampleDoseMatrix = [2 0 2; ...
                    1 2 0; ...
                    2 2 2; ...
                    4 4 4];
weights = [0.2; 0.3; 0.5];
analysis = matRad_expectedDoseDifferenceAnalysis(sampleDoseMatrix, referenceCube, ...
                                                 'scenWeights', weights);
fig = figure('Visible', 'off');
axesHandle = axes(fig);

matRad_plotExpectedDoseDifferenceAnalysis(analysis, struct(), {}, 1, ...
                                          'axesHandle', axesHandle);

imageHandle = findobj(axesHandle, 'Type', 'Image');
colorBarHandle = findobj(fig, 'Type', 'ColorBar');
assertEqual(numel(imageHandle), 1);
assertElementsAlmostEqual(get(imageHandle, 'CData'), analysis.signedExpectedDoseDifferenceCube, ...
                          'absolute', 1e-12);
assertElementsAlmostEqual(caxis(axesHandle), analysis.doseWindow, 'absolute', 1e-12);
assertEqual(get(get(colorBarHandle, 'YLabel'), 'String'), 'E[D - ref] [Gy]');

function test_plotUsesRbeDoseUnitLabel
figureCleaner = onCleanup(@() close('all'));
referenceCube = ones(2, 2);
sampleDoseMatrix = ones(4, 2);
analysis = matRad_expectedDoseDifferenceAnalysis(sampleDoseMatrix, referenceCube, ...
                                                 'quantity', 'RBExD');
fig = figure('Visible', 'off');
axesHandle = axes(fig);

matRad_plotExpectedDoseDifferenceAnalysis(analysis, struct(), {}, 1, ...
                                          'axesHandle', axesHandle);

colorBarHandle = findobj(fig, 'Type', 'ColorBar');
assertEqual(get(get(colorBarHandle, 'YLabel'), 'String'), 'E[D - ref] [Gy(RBE)]');

function test_invalidWeightsReturnSkippedAnalysis
referenceCube = ones(2, 2);
sampleDoseMatrix = ones(4, 3);

analysis = matRad_expectedDoseDifferenceAnalysis(sampleDoseMatrix, referenceCube, ...
                                                 'scenWeights', [1; 1]);

assertEqual(analysis.status, 'skippedInvalidInput');
assertFalse(isempty(strfind(analysis.reason, 'scenario weights')));
assertTrue(all(isnan(analysis.overReferenceProbabilityCube(:))));
