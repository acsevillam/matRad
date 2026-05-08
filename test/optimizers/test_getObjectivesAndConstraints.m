function test_suite = test_getObjectivesAndConstraints

test_functions = localfunctions();

initTestSuite;

function test_get_objectives_and_constraints_lists_prob2_functions

    classNames = matRad_getObjectivesAndConstraints();

    assertTrue(any(strcmp(classNames(1,:), ...
        'DoseObjectives.matRad_MeanVariance')));
    assertTrue(any(strcmp(classNames(1,:), ...
        'DoseConstraints.matRad_MinMaxMeanVariance')));
    assertTrue(any(strcmp(classNames(2,:), ...
        DoseObjectives.matRad_MeanVariance.name)));
    assertTrue(any(strcmp(classNames(2,:), ...
        DoseConstraints.matRad_MinMaxMeanVariance.name)));

function test_get_objectives_and_constraints_lists_interval_target_objective

    classNames = matRad_getObjectivesAndConstraints();

    assertTrue(any(strcmp(classNames(1,:), ...
        'DoseObjectives.matRad_SquaredBertoluzzaDeviation')));
    assertTrue(any(strcmp(classNames(2,:), ...
        DoseObjectives.matRad_SquaredBertoluzzaDeviation.name)));
