function test_suite = test_cheapCOWC

test_functions = localfunctions();

initTestSuite;

function test_cheapCowcDependencyIsAvailable
objective = DoseObjectives.matRad_SquaredDeviation(1, 0);
objective.robustness = 'c-COWC';

assertEqual(objective.robustness, 'c-COWC');
assertEqual(exist('gradest', 'file'), 2);
assertEqual(exist('derivest', 'file'), 2);

function test_cheapCowcDefaultRanksSelectWorstScenario
[optiProb, dij, cst, w] = helper_buildCheapCowcProblem([1 2 3], [0.2 0.3 0.5], 1, 1);

f = optiProb.matRad_objectiveFunction(w, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);

assertElementsAlmostEqual(f, 9, 'absolute', 1e-12);
assertElementsAlmostEqual(g, 18, 'absolute', 1e-6);

function test_cheapCowcRankWindowUsesProbabilityWeightedAverage
[optiProb, dij, cst, w] = helper_buildCheapCowcProblem([1 2 3], [0.2 0.3 0.5], 2, 3);

f = optiProb.matRad_objectiveFunction(w, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);

assertElementsAlmostEqual(f, 2.8, 'absolute', 1e-12);
assertElementsAlmostEqual(g, 5.6, 'absolute', 1e-6);

function test_cheapCowcSubsetsFullProbabilityVector
scenarioProb = [0.1 0.2 0.3 0.4];
[optiProb, dij, cst, w] = helper_buildCheapCowcProblem([0 2 0 3], scenarioProb, 2, 2);
optiProb.BP.scenarios = [2 4];

f = optiProb.matRad_objectiveFunction(w, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);

assertElementsAlmostEqual(f, 4, 'absolute', 1e-12);
assertElementsAlmostEqual(g, 8, 'absolute', 1e-6);

function test_cheapCowcRejectsInvalidRankWindow
[optiProb, dij, cst, w] = helper_buildCheapCowcProblem([1 2 3], [0.2 0.3 0.5], 2, 4);

assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w, dij, cst), ...
                      'matRad:Error');

function test_classicCowcKeepsWorstCaseBehavior
[optiProb, dij, cst, w] = helper_buildCowcProblem([1 2 3]);
optiProb.useMaxApprox = 'none';

f = optiProb.matRad_objectiveFunction(w, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);

assertElementsAlmostEqual(f, 9, 'absolute', 1e-12);
assertElementsAlmostEqual(g, 18, 'absolute', 1e-12);

function test_cheapCowcGradientMatchesCentralDifference
[optiProb, dij, cst, w] = helper_buildCheapCowcProblem([1 2 3], [0.2 0.3 0.5], 2, 3);

f = @(x) optiProb.matRad_objectiveFunction(x, dij, cst);
g = optiProb.matRad_objectiveGradient(w, dij, cst);
gNum = helper_centralDifference(f, w);

assertElementsAlmostEqual(g, gNum, 'absolute', 1e-5);

function [optiProb, dij, cst, w] = helper_buildCheapCowcProblem(doseValues, scenarioProb, p1, p2)
[optiProb, dij, cst, w] = helper_buildRobustProblem(doseValues, scenarioProb, 'c-COWC');
optiProb.p1 = p1;
optiProb.p2 = p2;

function [optiProb, dij, cst, w] = helper_buildCowcProblem(doseValues)
scenarioProb = ones(1, numel(doseValues)) ./ numel(doseValues);
[optiProb, dij, cst, w] = helper_buildRobustProblem(doseValues, scenarioProb, 'COWC');

function [optiProb, dij, cst, w] = helper_buildRobustProblem(doseValues, scenarioProb, robustness)
numScenarios = numel(doseValues);
w = 1;

backProjection = matRad_DoseProjection();
backProjection.scenarios = 1:numScenarios;
backProjection.scenarioProb = scenarioProb(:);
backProjection.nominalCtScenarios = 1;
optiProb = matRad_OptimizationProblem(backProjection);

dij.physicalDose = cell(1, numScenarios);
for i = 1:numScenarios
    dij.physicalDose{i} = sparse(doseValues(i));
end
dij.doseGrid.numOfVoxels = 1;
dij.totalNumOfBixels = 1;

objective = DoseObjectives.matRad_SquaredDeviation(1, 0);
objective.robustness = robustness;

cst = cell(1, 6);
cst{1, 3} = 'TARGET';
cst{1, 4} = {1};
cst{1, 5}.alphaX = [];
cst{1, 5}.betaX = [];
cst{1, 6} = {objective};

function g = helper_centralDifference(f, w)
delta = 1e-6;
g = zeros(size(w));

for i = 1:numel(w)
    step = zeros(size(w));
    step(i) = delta;
    g(i) = (f(w + step) - f(w - step)) / (2 * delta);
end
