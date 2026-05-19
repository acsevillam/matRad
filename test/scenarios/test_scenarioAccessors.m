function test_suite = test_scenarioAccessors

test_functions = localfunctions();

initTestSuite;

function test_scenarioAccessorsExposeRealizationMetadata
ct.numOfCtScen = 3;
rngState = rng;
randomCleaner = onCleanup(@() rng(rngState));
rng(0);
models = { ...
          matRad_NominalScenario(ct), ...
          matRad_WorstCaseScenarios(ct), ...
          matRad_ImportanceScenarios(ct), ...
          matRad_RandomScenarios(ct)};

for i = 1:numel(models)
    helper_assertScenarioAccessors(models{i});
end

function test_scenarioAccessorsResolveSparseCtScenarioIds
ct.numOfCtScen = 3;
model = matRad_NominalScenario(ct);
model.ctScenProb = [2 1];

assertEqual(model.scenarioIds(), 1);
assertEqual(model.numScenarios(), 1);
assertEqual(model.scenarioCtScenIds, 2);
assertEqual(model.getCtScenario(1), 2);
assertEqual(model.getDijScenarioIndex(1), 2);
assertEqual(model.getDijScenarioIndexBySubscripts(1, 1, 1, 'position'), 2);
assertEqual(model.getDijScenarioIndexBySubscripts(2, 1, 1, 'id'), 2);
assertTrue(model.isScenarioActiveBySubscripts(2, 1, 1, 'id'));
assertExceptionThrown(@() model.getDijScenarioIndexBySubscripts(2, 1, 1, ...
                                                                'position'), 'matRad:Error');

function test_scenarioAccessorsIdentifyNominalScenarios
ct.numOfCtScen = 3;
model = matRad_NominalScenario(ct);

assertEqual(model.getNominalScenarioIds(), model.scenarioIds());

function test_scenarioAccessorsExtractSparseCtScenario
ct.numOfCtScen = 3;
model = matRad_NominalScenario(ct);
model.ctScenProb = [2 1];

singleScenario = model.extractSingleScenario(1);

assertEqual(singleScenario.numScenarios(), 1);
assertEqual(singleScenario.getCtScenario(1), 2);
assertEqual(singleScenario.getDijScenarioIndex(1), 2);

function test_scenarioAccessorsDoNotPolluteStructCreation
ct.numOfCtScen = 3;
model = matRad_NominalScenario(ct);

modelStruct = struct('model', model.shortName, ...
                     'scenarioComponents', model.scenarioComponents, ...
                     'scenarioValueNames', {model.scenarioValueNames}, ...
                     'scenarioValues', model.scenarioValues, ...
                     'scenarioIdList', model.scenarioIdList, ...
                     'scenarioCtScenIds', model.scenarioCtScenIds, ...
                     'scenarioStorageSubscripts', model.scenarioStorageSubscripts);
assertTrue(isfield(modelStruct, 'scenarioComponents'));
createdModel = matRad_ScenarioModel.create(modelStruct, ct);

assertEqual(createdModel.numScenarios(), model.numScenarios());

function helper_assertScenarioAccessors(model)
ids = model.scenarioIds();
checkedRows = unique([1; model.totNumScen]);
expectedValueNames = {'setup.x', 'setup.y', 'setup.z', ...
                      'range.absolute', 'range.relative'};

assertEqual(ids, (1:model.totNumScen)');
assertEqual(model.numScenarios(), model.totNumScen);
assertEqual(model.scenarioCtScenIds, model.ctScenIx(:));
assertEqual(model.scenarioValueNames, expectedValueNames);
assertEqual({model.scenarioComponents.name}, expectedValueNames);
assertEqual(model.scenarioValues, model.scenForProb(:, 2:end));
assertEqual(model.getDijActiveMask(), model.scenMask);
assertEqual(model.getDijContainerSize(), size(model.scenMask));

for i = 1:numel(checkedRows)
    rowIx = checkedRows(i);
    scenarioId = ids(rowIx);
    scenario = model.getScenario(scenarioId);
    linearMask = model.linearMask(rowIx, :);

    assertEqual(scenario.id, scenarioId);
    assertEqual(scenario.ctScenId, model.ctScenIx(rowIx));
    assertEqual(scenario.values.setup_x, model.isoShift(rowIx, 1));
    assertEqual(scenario.values.setup_y, model.isoShift(rowIx, 2));
    assertEqual(scenario.values.setup_z, model.isoShift(rowIx, 3));
    assertEqual(scenario.values.range_absolute, model.absRangeShift(rowIx));
    assertEqual(scenario.values.range_relative, model.relRangeShift(rowIx));
    assertEqual(scenario.probability, model.scenProb(rowIx));
    assertEqual(scenario.weight, model.scenWeight(rowIx));
    assertEqual(model.getCtScenario(scenarioId), model.ctScenIx(rowIx));
    assertEqual(model.getSetupShift(scenarioId), model.isoShift(rowIx, :));
    assertEqual(model.getRangeShift(scenarioId), ...
                [model.absRangeShift(rowIx), model.relRangeShift(rowIx)]);
    assertEqual(model.getValue(scenarioId, 'setup.x'), model.isoShift(rowIx, 1));
    assertEqual(model.getValues(scenarioId, {'range.absolute', 'range.relative'}), ...
                [model.absRangeShift(rowIx), model.relRangeShift(rowIx)]);
    assertEqual(model.getDijScenarioIndex(scenarioId), ...
                model.sub2scenIx(linearMask(1), linearMask(2), linearMask(3)));
    assertEqual(model.getScenarioRowIndexFromDijIndex(model.getDijScenarioIndex(scenarioId)), ...
                rowIx);
    assertTrue(model.isScenarioActiveBySubscripts(linearMask(1), linearMask(2), ...
                                                  linearMask(3), 'id'));
end

inactiveScenIxs = find(~model.scenMask);
if ~isempty(inactiveScenIxs)
    assertTrue(isempty(model.getScenarioRowIndexFromDijIndex(inactiveScenIxs(1))));
end

assertExceptionThrown(@() model.getScenario(model.totNumScen + 1), ...
                      'matRad:Error');
assertExceptionThrown(@() model.getValue(ids(1), 'missing.component'), ...
                      'matRad:Error');
