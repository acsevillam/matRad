function test_suite = test_truncatedImportanceScenarios

test_functions=localfunctions();

initTestSuite;

function test_truncatedImportanceScenarioConstructor
    scenario = matRad_TruncatedImportanceScenarios();
    refScenario = matRad_ImportanceScenarios();

    assertTrue(isa(scenario,'matRad_TruncatedImportanceScenarios'));
    assertTrue(isa(scenario,'matRad_GriddedScenariosAbstract'));
    assertTrue(isa(scenario,'matRad_ScenarioModel'));
    assertEqual(scenario.name,'truncatedImpScen');

    % The truncated model starts from the native importance grid and removes
    % scenarios outside the configured uncertainty radius.
    assertTrue(scenario.totNumScen < refScenario.totNumScen);
    assertAllInsideTruncationRadius(scenario);
    assertElementsAlmostEqual(scenario.scenWeight,scenario.scenProb./sum(scenario.scenProb));

function test_truncatedImportanceScenarioConstructorWithCt
    n = 3;
    ct = struct('numOfCtScen',n);
    scenario = matRad_TruncatedImportanceScenarios(ct);

    assertEqual(scenario.ctScenProb,[(1:n)' ones(n,1)./n]);
    assertEqual(scenario.numOfCtScen,n);
    assertEqual(size(scenario.scenMask,1),n);
    assertEqual(size(scenario.linearMask,1),scenario.totNumScen);

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
    assertAllInsideTruncationRadius(scenario);

    [tmp(:,1),tmp(:,2),tmp(:,3)] = ind2sub(size(scenario.scenMask),find(scenario.scenMask));
    assertEqual(tmp,scenario.linearMask);

function test_truncatedImportanceScenarioFactory
    model = matRad_multScen([],'truncatedImpScen');
    assertTrue(isa(model,'matRad_TruncatedImportanceScenarios'));
    assertEqual(model.name,'truncatedImpScen');

function assertAllInsideTruncationRadius(scenario)
    scenValues = scenario.scenForProb(:,2:6);
    scenScale = [scenario.shiftSD scenario.rangeAbsSD scenario.rangeRelSD./100];
    scenScale(scenScale == 0) = eps;
    normalizedRadius = sqrt(sum(bsxfun(@rdivide,scenValues,scenScale).^2,2));
    tolerance = 100*eps(max(1,scenario.wcSigma));

    assertTrue(all(normalizedRadius <= scenario.wcSigma + tolerance));
