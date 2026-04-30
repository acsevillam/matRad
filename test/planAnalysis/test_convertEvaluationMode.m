function test_suite = test_convertEvaluationMode

test_functions = localfunctions();

initTestSuite;

function test_convert_to_total
    pln.numOfFractions = 39;

    [value,evaluationMode,evaluationScale] = matRad_convertToEvaluationMode( ...
        [0 2],pln,'total');

    assertElementsAlmostEqual(value,[0 78]);
    assertEqual(evaluationMode,'total');
    assertEqual(evaluationScale,39);

function test_convert_from_total
    pln.numOfFractions = 39;

    [value,evaluationMode,evaluationScale] = matRad_convertFromEvaluationMode( ...
        [0 78],pln,'total');

    assertElementsAlmostEqual(value,[0 2]);
    assertEqual(evaluationMode,'total');
    assertEqual(evaluationScale,39);

function test_convert_from_per_fraction
    pln.numOfFractions = 39;

    [value,evaluationMode,evaluationScale] = matRad_convertFromEvaluationMode( ...
        [0 2],pln,'perFraction');

    assertElementsAlmostEqual(value,[0 2]);
    assertEqual(evaluationMode,'perFraction');
    assertEqual(evaluationScale,1);

function test_convert_from_empty_value
    pln.numOfFractions = 39;

    value = matRad_convertFromEvaluationMode([],pln,'total');

    assertTrue(isempty(value));
