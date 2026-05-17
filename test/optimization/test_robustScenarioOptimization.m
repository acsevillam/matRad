function test_suite = test_robustScenarioOptimization

test_functions = localfunctions();

initTestSuite;

function test_prob2MeanVarianceValueAndGradient
[optiProb, dij, cst, w] = helper_prob2MeanVarianceProblem();

f = optiProb.matRad_objectiveFunction(w, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);

assertElementsAlmostEqual(f, 12, 'absolute', 1e-12);
assertElementsAlmostEqual(g, [9; 7.5], 'absolute', 1e-12);

function test_prob2MeanVarianceGradientMatchesCentralDifference
[optiProb, dij, cst, w] = helper_prob2MeanVarianceProblem();

f = @(x) optiProb.matRad_objectiveFunction(x, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);
gNum = helper_centralDifference(f, w);

assertElementsAlmostEqual(g, gNum, 'absolute', 1e-6);

function test_prob2StandardObjectiveUsesExpectedDose
[optiProb, dij, cst, w] = helper_prob2ExpectedDoseProblem();

f = optiProb.matRad_objectiveFunction(w, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);

assertElementsAlmostEqual(f, 9, 'absolute', 1e-12);
assertElementsAlmostEqual(g, [0; 12], 'absolute', 1e-12);

function test_probObjectiveUsesExpectedDoseAndOmega
[optiProb, dij, cst, w] = helper_probExpectedDoseProblem();

f = optiProb.matRad_objectiveFunction(w, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);

assertElementsAlmostEqual(f, 17, 'absolute', 1e-12);
assertElementsAlmostEqual(g, [6; 17], 'absolute', 1e-12);

function test_probGradientMatchesCentralDifference
[optiProb, dij, cst, w] = helper_probExpectedDoseProblem();

f = @(x) optiProb.matRad_objectiveFunction(x, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);
gNum = helper_centralDifference(f, w);

assertElementsAlmostEqual(g, gNum, 'absolute', 1e-6);

function test_probConstraintUsesExpectedDose
[optiProb, dij, cst, w] = helper_probExpectedDoseConstraintProblem();

c = optiProb.matRad_constraintFunctions(w, dij, cst);
jacob = optiProb.matRad_constraintJacobian(w, dij, cst);
jacobStruct = optiProb.matRad_getJacobianStructure(w, dij, cst);

assertElementsAlmostEqual(c, 2.5, 'absolute', 1e-12);
assertElementsAlmostEqual(full(jacob), [0.5 1], 'absolute', 1e-12);
assertEqual(full(jacobStruct), [1 1]);

function test_prob2MeanVarianceConstraintValueAndJacobian
[optiProb, dij, cst, w] = helper_prob2MeanVarianceConstraintProblem();

c = optiProb.matRad_constraintFunctions(w, dij, cst);
jacob = optiProb.matRad_constraintJacobian(w, dij, cst);
jacobStruct = optiProb.matRad_getJacobianStructure(w, dij, cst);

assertElementsAlmostEqual(c, 4, 'absolute', 1e-12);
assertElementsAlmostEqual(full(jacob), [3 2.5], 'absolute', 1e-12);
assertEqual(full(jacobStruct), [1 1]);

function test_prob2ConstraintsUseReferenceScenarioContour
[optiProb, dij, cst, w] = helper_prob2MeanVarianceConstraintProblem();
optiProb.dij_prob.refScen = 2;
cst{1, 4} = {[], [1; 2]};

c = optiProb.matRad_constraintFunctions(w, dij, cst);
jacob = optiProb.matRad_constraintJacobian(w, dij, cst);
jacobStruct = optiProb.matRad_getJacobianStructure(w, dij, cst);

assertElementsAlmostEqual(c, 4, 'absolute', 1e-12);
assertElementsAlmostEqual(full(jacob), [3 2.5], 'absolute', 1e-12);
assertEqual(full(jacobStruct), [1 1]);

function test_mixedConstraintJacobianPreservesConstraintOrder
[optiProb, dij, cst, w] = helper_mixedConstraintProblem(false);

c = optiProb.matRad_constraintFunctions(w, dij, cst);
jacob = optiProb.matRad_constraintJacobian(w, dij, cst);
jacobStruct = optiProb.matRad_getJacobianStructure(w, dij, cst);

assertElementsAlmostEqual(c, [2.5; 4], 'absolute', 1e-12);
assertElementsAlmostEqual(full(jacob), [0.5 1; 3 2.5], ...
                          'absolute', 1e-12);
assertEqual(full(jacobStruct), [1 1; 1 1]);

[optiProb, dij, cst, w] = helper_mixedConstraintProblem(true);

c = optiProb.matRad_constraintFunctions(w, dij, cst);
jacob = optiProb.matRad_constraintJacobian(w, dij, cst);
jacobStruct = optiProb.matRad_getJacobianStructure(w, dij, cst);

assertElementsAlmostEqual(c, [4; 2.5], 'absolute', 1e-12);
assertElementsAlmostEqual(full(jacob), [3 2.5; 0.5 1], ...
                          'absolute', 1e-12);
assertEqual(full(jacobStruct), [1 1; 1 1]);

function test_prob2RejectsMissingPayload
[optiProb, dij, cst, w] = helper_prob2MeanVarianceProblem();
optiProb.dij_prob = struct();

assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w, dij, cst), ...
                      'matRad:Error');

function test_probRejectsMissingPayload
[optiProb, dij, cst, w] = helper_probExpectedDoseProblem();
optiProb.dij_prob = struct();

assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w, dij, cst), ...
                      'matRad:Error');

function test_prob2RejectsQuantityMismatch
[optiProb, dij, cst, w] = helper_prob2ExpectedDoseProblem();
optiProb.quantityOpt = 'RBExDose';
optiProb.dij_prob.quantity = 'physicalDose';

assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w, dij, cst), ...
                      'matRad:Error');

function test_probRejectsMeanVarianceConstraint
assertExceptionThrown(@() helper_assignProbRobustnessToMeanVarianceConstraint(), ...
                      'matRad:Error');

function test_prob2RejectsMissingOmegaForMeanVarianceConstraint
[optiProb, dij, cst, w] = helper_prob2MeanVarianceConstraintProblem();
optiProb.dij_prob.Omega = {};

assertExceptionThrown(@() optiProb.matRad_constraintFunctions(w, dij, cst), ...
                      'matRad:Error');

function test_intervalTargetBertoluzzaValueAndGradient
[optiProb, dij, cst, w] = helper_intervalTargetProblem('INTERVAL2');

f = optiProb.matRad_objectiveFunction(w, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);

assertElementsAlmostEqual(f, 8.1, 'absolute', 1e-12);
assertElementsAlmostEqual(g, [2.4; 11.4], 'absolute', 1e-10);

function test_intervalTargetGradientMatchesCentralDifference
[optiProb, dij, cst, w] = helper_intervalTargetProblem('INTERVAL3');

f = @(x) optiProb.matRad_objectiveFunction(x, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);
gNum = helper_centralDifference(f, w);

assertElementsAlmostEqual(g, gNum, 'absolute', 1e-6);

function test_interval3OarUsesCenterPlusRadius
[optiProb, dij, cst, w, expectedDose] = helper_interval3OARProblem();
objective = cst{1, 6}{1};

f = optiProb.matRad_objectiveFunction(w, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);
gNum = helper_centralDifference(@(x) optiProb.matRad_objectiveFunction(x, dij, cst), w);
expectedF = objective.penalty * objective.computeDoseObjectiveFunction(expectedDose);

assertElementsAlmostEqual(f, expectedF, 'absolute', 1e-12);
assertElementsAlmostEqual(g, gNum, 'absolute', 1e-6);

function test_intervalRejectsMissingPayload
[optiProb, dij, cst, w] = helper_intervalTargetProblem('INTERVAL2');
optiProb.dij_interval = struct();

assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w, dij, cst), ...
                      'matRad:Error');

function test_intervalRejectsQuantityMismatch
[optiProb, dij, cst, w] = helper_intervalTargetProblem('INTERVAL2');
optiProb.quantityOpt = 'RBExDose';
optiProb.dij_interval.quantity = 'physicalDose';

assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w, dij, cst), ...
                      'matRad:Error');

function test_nominalObjectiveRemainsUnchanged
[optiProb, dij, cst, w] = helper_nominalProblem();

f = optiProb.matRad_objectiveFunction(w, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);

assertElementsAlmostEqual(f, 8, 'absolute', 1e-12);
assertElementsAlmostEqual(g, [8; 8], 'absolute', 1e-12);

function [optiProb, dij, cst, w] = helper_prob2MeanVarianceProblem()
[optiProb, dij, cst, w] = helper_baseProb2Problem();
objective = DoseObjectives.matRad_MeanVariance(3);
objective.robustness = 'PROB2';
cst{1, 6} = {objective};

function [optiProb, dij, cst, w] = helper_prob2ExpectedDoseProblem()
[optiProb, dij, cst, w] = helper_baseProb2Problem();
objective = DoseObjectives.matRad_SquaredDeviation(2, 1);
objective.robustness = 'PROB2';
cst{1, 6} = {objective};

function [optiProb, dij, cst, w] = helper_probExpectedDoseProblem()
[optiProb, dij, cst, w] = helper_baseProb2Problem();
objective = DoseObjectives.matRad_SquaredDeviation(2, 1);
objective.robustness = 'PROB';
cst{1, 6} = {objective};

function [optiProb, dij, cst, w] = helper_probExpectedDoseConstraintProblem()
[optiProb, dij, cst, w] = helper_baseProb2Problem();
constraint = DoseConstraints.matRad_MinMaxMeanDose(0, 10);
constraint.robustness = 'PROB';
cst{1, 6} = {constraint};

function [optiProb, dij, cst, w] = helper_prob2MeanVarianceConstraintProblem()
[optiProb, dij, cst, w] = helper_baseProb2Problem();
constraint = DoseConstraints.matRad_MinMaxMeanVariance(0, 10);
constraint.robustness = 'PROB2';
cst{1, 6} = {constraint};

function [optiProb, dij, cst, w] = helper_mixedConstraintProblem(prob2First)
[optiProb, dij, cst, w] = helper_baseProb2Problem();
probConstraint = DoseConstraints.matRad_MinMaxMeanDose(0, 10);
probConstraint.robustness = 'PROB';
prob2Constraint = DoseConstraints.matRad_MinMaxMeanVariance(0, 10);
prob2Constraint.robustness = 'PROB2';

if prob2First
    cst{1, 6} = {prob2Constraint, probConstraint};
else
    cst{1, 6} = {probConstraint, prob2Constraint};
end

function helper_assignProbRobustnessToMeanVarianceConstraint()
constraint = DoseConstraints.matRad_MinMaxMeanVariance(0, 10);
constraint.robustness = 'PROB';

function [optiProb, dij, cst, w] = helper_baseProb2Problem()
w = [1; 2];
expected = sparse([1 0; 0 2; 3 0]);
Omega = sparse([2 0.5; 0.5 1]);

optiProb = helper_createOptimizationProblem();
optiProb.quantityOpt = 'physicalDose';
optiProb.dij_prob.expected = expected;
optiProb.dij_prob.Omega = cell(1, 1);
optiProb.dij_prob.Omega{1} = Omega;
optiProb.dij_prob.voiSubIx = cell(1, 1);
optiProb.dij_prob.voiSubIx{1} = [1; 2];
optiProb.dij_prob.quantity = 'physicalDose';
optiProb.dij_prob.refScen = 1;

dij = helper_createDij(expected);
cst = helper_createCst('TARGET', [1; 2]);

function [optiProb, dij, cst, w] = helper_intervalTargetProblem(robustness)
w = [1; 2];
optiProb = helper_createOptimizationProblem();
optiProb.quantityOpt = 'physicalDose';
optiProb.theta1 = 0.4;
optiProb.dij_interval.center = sparse([1 0; 0 2; 3 0]);
optiProb.dij_interval.radius = sparse([2 0.5; 0.5 1]);
optiProb.dij_interval.quantity = 'physicalDose';
optiProb.dij_interval.refScen = 1;

dij = helper_createDij(sparse(3, 2));
objective = DoseObjectives.matRad_SquaredBertoluzzaDeviation(3, 1);
objective.robustness = robustness;
cst = helper_createCst('TARGET', [1; 2]);
cst{1, 6} = {objective};

function [optiProb, dij, cst, w, expectedDose] = helper_interval3OARProblem()
w = [3; 4];
optiProb = helper_createOptimizationProblem();
optiProb.quantityOpt = 'physicalDose';
optiProb.theta2 = 0.5;
optiProb.dij_interval.center = sparse([1 0; 0 1; 0 0]);
optiProb.dij_interval.radius = sparse(2, 2);
optiProb.dij_interval.OARSubIx = [1; 2];
optiProb.dij_interval.OARRadiusFactor = {sparse(diag([2 1])); ...
                                         sparse(diag([1 3]))};
optiProb.dij_interval.OARRadiusRank = [2; 2];
optiProb.dij_interval.quantity = 'physicalDose';
optiProb.dij_interval.refScen = 1;

dCenter = optiProb.dij_interval.center([1; 2], :) * w;
dRadius = [norm(optiProb.dij_interval.OARRadiusFactor{1}' * w); ...
           norm(optiProb.dij_interval.OARRadiusFactor{2}' * w)];
expectedDose = dCenter + optiProb.theta2 * dRadius;

dij = helper_createDij(sparse(3, 2));
objective = DoseObjectives.matRad_SquaredDeviation(2, 1);
objective.robustness = 'INTERVAL3';
cst = helper_createCst('OAR', [1; 2]);
cst{1, 6} = {objective};

function [optiProb, dij, cst, w] = helper_nominalProblem()
w = [1; 2];
doseMatrix = sparse([1 1; 1 1]);
optiProb = helper_createOptimizationProblem();
dij = helper_createDij(doseMatrix);
objective = DoseObjectives.matRad_SquaredDeviation(2, 1);
cst = helper_createCst('TARGET', [1; 2]);
cst{1, 6} = {objective};

function optiProb = helper_createOptimizationProblem()
backProjection = matRad_DoseProjection();
backProjection.scenarios = 1;
backProjection.scenarioProb = 1;
backProjection.nominalCtScenarios = 1;
optiProb = matRad_OptimizationProblem(backProjection);

function dij = helper_createDij(doseMatrix)
dij.physicalDose = {doseMatrix};
dij.doseGrid.numOfVoxels = size(doseMatrix, 1);
dij.totalNumOfBixels = size(doseMatrix, 2);

function cst = helper_createCst(voiType, voxels)
cst = cell(1, 6);
cst{1, 3} = voiType;
cst{1, 4} = {voxels(:)};
cst{1, 5}.alphaX = [];
cst{1, 5}.betaX = [];
cst{1, 6} = {};

function g = helper_centralDifference(f, w)
delta = 1e-6;
g = zeros(size(w));

for i = 1:numel(w)
    step = zeros(size(w));
    step(i) = delta;
    g(i) = (f(w + step) - f(w - step)) / (2 * delta);
end
