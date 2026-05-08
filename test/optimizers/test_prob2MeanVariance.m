function test_suite = test_prob2MeanVariance

test_functions = localfunctions();

initTestSuite;

function test_mean_variance_objective_value_and_gradient
    [optiProb,dij,cst,w] = buildMeanVarianceProblem();

    f = optiProb.matRad_objectiveFunction(w,dij,cst);
    g = optiProb.matRad_objectiveGradient(w,dij,cst);

    assertElementsAlmostEqual(f,12,'absolute',1e-12);
    assertElementsAlmostEqual(g,[9;7.5],'absolute',1e-12);

function test_mean_variance_gradient_matches_central_difference
    [optiProb,dij,cst,w] = buildMeanVarianceProblem();

    f = @(x) optiProb.matRad_objectiveFunction(x,dij,cst);
    g = optiProb.matRad_objectiveGradient(w,dij,cst);
    gNum = centralDifference(f,w);

    assertElementsAlmostEqual(g,gNum,'absolute',1e-6);

function test_standard_prob2_objective_uses_expected_dose
    [optiProb,dij,cst,w] = buildStandardDoseProblem();

    f = optiProb.matRad_objectiveFunction(w,dij,cst);
    g = optiProb.matRad_objectiveGradient(w,dij,cst);

    assertElementsAlmostEqual(f,9,'absolute',1e-12);
    assertElementsAlmostEqual(g,[0;12],'absolute',1e-12);

function test_mean_variance_constraint_value_and_jacobian
    [optiProb,dij,cst,w] = buildMeanVarianceConstraintProblem();

    c = optiProb.matRad_constraintFunctions(w,dij,cst);
    jacob = optiProb.matRad_constraintJacobian(w,dij,cst);

    assertElementsAlmostEqual(c,4,'absolute',1e-12);
    assertElementsAlmostEqual(full(jacob),[3 2.5],'absolute',1e-12);

function test_prob2_constraints_use_reference_scenario_contour
    [optiProb,dij,cst,w] = buildMeanVarianceConstraintProblem();
    optiProb.dij_prob2.refScen = 2;
    cst{1,4} = {[],[1;2]};

    c = optiProb.matRad_constraintFunctions(w,dij,cst);
    jacob = optiProb.matRad_constraintJacobian(w,dij,cst);
    jacobStruct = optiProb.matRad_getJacobianStructure(zeros(size(w)),dij,cst);

    assertElementsAlmostEqual(c,4,'absolute',1e-12);
    assertElementsAlmostEqual(full(jacob),[3 2.5],'absolute',1e-12);
    assertEqual(full(jacobStruct),[1 1]);

function test_mixed_constraint_jacobian_preserves_constraint_order
    [optiProb,dij,cst,w] = buildMixedConstraintProblem(false);

    c = optiProb.matRad_constraintFunctions(w,dij,cst);
    jacob = optiProb.matRad_constraintJacobian(w,dij,cst);
    jacobStruct = optiProb.matRad_getJacobianStructure(zeros(size(w)),dij,cst);

    assertElementsAlmostEqual(c,[2.5;4],'absolute',1e-12);
    assertElementsAlmostEqual(full(jacob),[0.5 1; 3 2.5], ...
        'absolute',1e-12);
    assertEqual(full(jacobStruct),[1 1; 1 1]);

    [optiProb,dij,cst,w] = buildMixedConstraintProblem(true);

    c = optiProb.matRad_constraintFunctions(w,dij,cst);
    jacob = optiProb.matRad_constraintJacobian(w,dij,cst);
    jacobStruct = optiProb.matRad_getJacobianStructure(zeros(size(w)),dij,cst);

    assertElementsAlmostEqual(c,[4;2.5],'absolute',1e-12);
    assertElementsAlmostEqual(full(jacob),[3 2.5; 0.5 1], ...
        'absolute',1e-12);
    assertEqual(full(jacobStruct),[1 1; 1 1]);

function test_prob2_rejects_missing_dij_prob2
    [optiProb,dij,cst,w] = buildMeanVarianceProblem();
    optiProb.dij_prob2 = struct();

    assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w,dij,cst), ...
        'matRad:Error');

function test_prob2_public_result_stats
    [optiProb,dij,cst,w] = buildMeanVarianceProblem();

    stats = optiProb.GetResultProbabilistic(w,dij,cst,1);

    assertElementsAlmostEqual(stats.dExp,[1;4],'absolute',1e-12);
    assertElementsAlmostEqual(stats.omegaW,[3;2.5],'absolute',1e-12);
    assertElementsAlmostEqual(stats.meanVariance,4,'absolute',1e-12);
    assertElementsAlmostEqual(stats.gradMeanVariance,[3;2.5],'absolute',1e-12);

function [optiProb,dij,cst,w] = buildMeanVarianceProblem()
    [optiProb,dij,cst,w] = buildBaseProb2Problem();
    objective = DoseObjectives.matRad_MeanVariance(3);
    objective.robustness = 'PROB2';
    cst{1,6} = {objective};

function [optiProb,dij,cst,w] = buildStandardDoseProblem()
    [optiProb,dij,cst,w] = buildBaseProb2Problem();
    objective = DoseObjectives.matRad_SquaredDeviation(2,1);
    objective.robustness = 'PROB2';
    cst{1,6} = {objective};

function [optiProb,dij,cst,w] = buildMeanVarianceConstraintProblem()
    [optiProb,dij,cst,w] = buildBaseProb2Problem();
    constraint = DoseConstraints.matRad_MinMaxMeanVariance(0,10);
    constraint.robustness = 'PROB2';
    cst{1,6} = {constraint};

function [optiProb,dij,cst,w] = buildMixedConstraintProblem(prob2First)
    [optiProb,dij,cst,w] = buildBaseProb2Problem();

    nominalConstraint = DoseConstraints.matRad_MinMaxMeanDose(0,10);
    prob2Constraint = DoseConstraints.matRad_MinMaxMeanVariance(0,10);
    prob2Constraint.robustness = 'PROB2';

    if prob2First
        cst{1,6} = {prob2Constraint,nominalConstraint};
    else
        cst{1,6} = {nominalConstraint,prob2Constraint};
    end

function [optiProb,dij,cst,w] = buildBaseProb2Problem()
    w = [1;2];
    expected = sparse([1 0; 0 2; 3 0]);
    Omega = sparse([2 0.5; 0.5 1]);

    backProjection = matRad_DoseProjection();
    backProjection.scenarios = 1;
    backProjection.scenarioProb = 1;
    backProjection.nominalCtScenarios = 1;

    optiProb = matRad_OptimizationProblem(backProjection);
    optiProb.quantityOpt = 'physicalDose';
    optiProb.dij_prob2.expected = expected;
    optiProb.dij_prob2.Omega = cell(1,1);
    optiProb.dij_prob2.Omega{1} = Omega;
    optiProb.dij_prob2.voiSubIx = cell(1,1);
    optiProb.dij_prob2.voiSubIx{1} = [1;2];
    optiProb.dij_prob2.quantity = 'physicalDose';
    optiProb.dij_prob2.refScen = 1;

    dij.physicalDose = {expected};
    dij.doseGrid.numOfVoxels = 3;
    dij.totalNumOfBixels = 2;

    cst = cell(1,6);
    cst{1,3} = 'TARGET';
    cst{1,4} = {[1;2]};
    cst{1,5}.alphaX = [];
    cst{1,5}.betaX = [];
    cst{1,6} = {};

function g = centralDifference(f,w)
    delta = 1e-6;
    g = zeros(size(w));

    for i = 1:numel(w)
        step = zeros(size(w));
        step(i) = delta;
        g(i) = (f(w + step) - f(w - step)) / (2*delta);
    end
