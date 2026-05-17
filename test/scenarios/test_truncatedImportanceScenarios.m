function test_suite = test_truncatedImportanceScenarios

test_functions = localfunctions();

initTestSuite;

function test_truncatedImportanceScenarioConstructor
scenario = matRad_TruncatedImportanceScenarios();
refScenario = matRad_ImportanceScenarios();

assertTrue(isa(scenario, 'matRad_TruncatedImportanceScenarios'));
assertTrue(isa(scenario, 'matRad_GriddedScenariosAbstract'));
assertTrue(isa(scenario, 'matRad_ScenarioModel'));
assertEqual(scenario.shortName, 'truncatedImpScen');
assertTrue(scenario.totNumScen < refScenario.totNumScen);
helper_assertAllInsideTruncationRadius(scenario);
expectedWeight = scenario.scenProb ./ sum(scenario.scenProb);
assertElementsAlmostEqual(scenario.scenWeight, expectedWeight);

function test_truncatedImportanceScenarioConstructorWithCt
n = 3;
ct = struct('numOfCtScen', n);
scenario = matRad_TruncatedImportanceScenarios(ct);

assertEqual(scenario.ctScenProb, [(1:n)' ones(n, 1) ./ n]);
assertEqual(scenario.numOfCtScen, n);
assertEqual(size(scenario.scenMask, 1), n);
assertEqual(size(scenario.linearMask, 1), scenario.totNumScen);

function test_truncatedImportanceScenarioRemovesCombinedCorners
scenario = matRad_TruncatedImportanceScenarios();
refScenario = matRad_ImportanceScenarios();

scenario.numOfSetupGridPoints = 3;
scenario.numOfRangeGridPoints = 3;
scenario.combineRange = false;
scenario.combinations = 'all';

refScenario.numOfSetupGridPoints = scenario.numOfSetupGridPoints;
refScenario.numOfRangeGridPoints = scenario.numOfRangeGridPoints;
refScenario.combineRange = scenario.combineRange;
refScenario.combinations = scenario.combinations;

assertTrue(scenario.totNumScen < refScenario.totNumScen);
helper_assertAllInsideTruncationRadius(scenario);

[tmp(:, 1), tmp(:, 2), tmp(:, 3)] = ind2sub(size(scenario.scenMask), find(scenario.scenMask));
assertEqual(tmp, scenario.linearMask);

function test_truncatedImportanceScenarioFactory
model = matRad_createScenarioModel([], 'truncatedImpScen');

assertTrue(isa(model, 'matRad_TruncatedImportanceScenarios'));
assertEqual(model.shortName, 'truncatedImpScen');

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
