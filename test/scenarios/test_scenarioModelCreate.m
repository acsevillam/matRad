function test_suite = test_scenarioModelCreate

test_functions = localfunctions();

initTestSuite;

function test_scenarioModelCreateWithCt
ct.numOfCtScen = 1;

helper_assertFactoryCreatesModel(ct, 'nomScen', 'matRad_NominalScenario');
helper_assertFactoryCreatesModel(ct, 'wcScen', 'matRad_WorstCaseScenarios');
helper_assertFactoryCreatesModel(ct, 'rndScen', 'matRad_RandomScenarios');
helper_assertFactoryCreatesModel(ct, 'impScen', 'matRad_ImportanceScenarios');

function test_scenarioModelCreateWithEmptyCt
ct = [];

helper_assertFactoryCreatesModel(ct, 'nomScen', 'matRad_NominalScenario');
helper_assertFactoryCreatesModel(ct, 'wcScen', 'matRad_WorstCaseScenarios');
helper_assertFactoryCreatesModel(ct, 'rndScen', 'matRad_RandomScenarios');
helper_assertFactoryCreatesModel(ct, 'impScen', 'matRad_ImportanceScenarios');

function test_scenarioModelCreateReturnsExistingInstance
model = matRad_NominalScenario();

createdModel = matRad_ScenarioModel.create(model);

assertTrue(createdModel == model);

function helper_assertFactoryCreatesModel(ct, modelName, className)

model = matRad_ScenarioModel.create(modelName, ct);

assertTrue(isa(model, className));
