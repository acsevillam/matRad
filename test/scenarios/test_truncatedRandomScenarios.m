function test_suite = test_truncatedRandomScenarios

test_functions = localfunctions();

initTestSuite;

function test_truncatedRandomScenarioConstructor
scenario = matRad_TruncatedRandomScenarios();

assertTrue(isa(scenario, 'matRad_TruncatedRandomScenarios'));
assertTrue(isa(scenario, 'matRad_RandomScenariosAbstract'));
assertTrue(isa(scenario, 'matRad_ScenarioModel'));
assertEqual(scenario.shortName, 'truncatedRndScen');
assertEqual(scenario.numScenarios(), scenario.nSamples);
helper_assertAllInsideTruncationRadius(scenario);

function test_truncatedRandomScenarioConstructorWithCt
n = 3;
ct = struct('numOfCtScen', n);
scenario = matRad_TruncatedRandomScenarios(ct);

assertEqual(scenario.ctScenProb, [(1:n)' ones(n, 1) ./ n]);
assertEqual(scenario.numOfCtScen, n);
assertEqual(scenario.numScenarios(), scenario.nSamples * n);
assertEqual(size(scenario.scenMask, 1), n);
helper_assertAllInsideTruncationRadius(scenario);

function test_truncatedRandomScenarioFactory
model = matRad_createScenarioModel([], 'truncatedRndScen');

assertTrue(isa(model, 'matRad_TruncatedRandomScenarios'));
assertEqual(model.shortName, 'truncatedRndScen');

function test_truncatedRandomScenarioReproducibleWithSeed
scenarioA = matRad_TruncatedRandomScenarios();
scenarioB = matRad_TruncatedRandomScenarios();

scenarioA.randomSeed = 42;
scenarioB.randomSeed = 42;
scenarioA.nSamples = 20;
scenarioB.nSamples = 20;

assertElementsAlmostEqual(scenarioA.scenarioValues, scenarioB.scenarioValues);
assertElementsAlmostEqual(scenarioA.scenWeight, scenarioB.scenWeight);

function test_truncatedRandomScenarioSupportsCompactAngularStorage
scenario = matRad_TruncatedRandomScenarios();

scenario.numOfBeams = 2;
scenario.scenarioDimensionActive = {'ct', 'setup', 'range', 'gantry'};
scenario.randomSeed = 7;
scenario.nSamples = 20;

assertEqual(scenario.scenarioStoragePolicy, 'compact-realization');
assertEqual(size(scenario.scenarioStorageSubscripts, 2), 2);
assertTrue(any(abs(scenario.gantryAngleOffset(:)) > 0));
assertEqual(size(scenario.gantryAngleOffset, 2), scenario.numOfBeams);
helper_assertAllInsideTruncationRadius(scenario);

function test_truncatedRandomScenarioZeroRadiusIsNominal
scenario = matRad_TruncatedRandomScenarios();

scenario.wcSigma = 0;

assertTrue(all(abs(scenario.scenarioValues(:)) <= eps));
assertTrue(all(abs(scenario.isoShift(:)) <= eps));
assertTrue(all(abs(scenario.absRangeShift(:)) <= eps));
assertTrue(all(abs(scenario.relRangeShift(:)) <= eps));

function helper_assertAllInsideTruncationRadius(scenario)
scenarioValues = scenario.scenarioValues;
activeIx = [scenario.scenarioComponents.active];
normalizedRadius = zeros(size(scenarioValues, 1), 1);

if any(activeIx)
    scenarioScale = [scenario.scenarioComponents(activeIx).scale];
    scaledValues = bsxfun(@rdivide, scenarioValues(:, activeIx), scenarioScale);
    normalizedRadius = sqrt(sum(scaledValues.^2, 2));
end

tolerance = 100 * eps(max(1, scenario.wcSigma));

assertTrue(all(normalizedRadius <= scenario.wcSigma + tolerance));
