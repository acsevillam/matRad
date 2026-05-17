function test_suite = test_convertEvaluationMode

test_functions = localfunctions();

initTestSuite;

function test_convertToTotal
pln.numOfFractions = 39;

[value, evaluationMode, evaluationScale] = matRad_convertToEvaluationMode( ...
                                                                          [0 2], pln, 'total');

assertElementsAlmostEqual(value, [0 78]);
assertEqual(evaluationMode, 'total');
assertEqual(evaluationScale, 39);

function test_convertFromTotal
pln.numOfFractions = 39;

[value, evaluationMode, evaluationScale] = matRad_convertFromEvaluationMode( ...
                                                                            [0 78], pln, 'total');

assertElementsAlmostEqual(value, [0 2]);
assertEqual(evaluationMode, 'total');
assertEqual(evaluationScale, 39);

function test_convertFromPerFraction
pln.numOfFractions = 39;

[value, evaluationMode, evaluationScale] = matRad_convertFromEvaluationMode( ...
                                                                            [0 2], pln, 'perFraction');

assertElementsAlmostEqual(value, [0 2]);
assertEqual(evaluationMode, 'perFraction');
assertEqual(evaluationScale, 1);

function test_convertFromEmptyValue
pln.numOfFractions = 39;

value = matRad_convertFromEvaluationMode([], pln, 'total');

assertTrue(isempty(value));

function test_invalidEvaluationModeFails
pln.numOfFractions = 39;

assertExceptionThrown(@() matRad_convertToEvaluationMode(1, pln, 'invalid'), ...
                      'matRad:Error');
