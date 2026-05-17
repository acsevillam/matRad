function test_suite = test_scenarioDimensionSupport

test_functions = localfunctions();

initTestSuite;

function test_pencilBeamSupportsAngularRandomScenariosInCompactStorage
rngState = rng;
rngCleaner = onCleanup(@() rng(rngState));
rng(0);

testData = load('protons_testData.mat');
testData.pln.propDoseCalc.engine = 'HongPB';
testData.pln.multScen = helper_createAngularRandomScenario(testData.ct, testData.stf);

dij = matRad_calcDoseInfluence(testData.ct, testData.cst, testData.stf, testData.pln);

multScen = testData.pln.multScen;
scenarioIds = multScen.scenarioIds();
activeMask = multScen.getDijActiveMask();

assertEqual(multScen.scenarioStoragePolicy, 'compact-realization');
assertEqual(multScen.getDijContainerSize(), [multScen.numScenarios() 1]);
assertEqual(numel(find(activeMask)), multScen.numScenarios());
assertEqual(dij.numOfScenarios, multScen.numScenarios());
assertEqual(numel(dij.physicalDose), numel(activeMask));

for i = 1:numel(scenarioIds)
    dijScenarioIx = multScen.getDijScenarioIndex(scenarioIds(i));
    assertTrue(activeMask(dijScenarioIx));
    assertFalse(isempty(dij.physicalDose{dijScenarioIx}));
    assertEqual(size(dij.physicalDose{dijScenarioIx}, 2), dij.totalNumOfBixels);
end

function test_nonMigratedEngineRejectsAngularScenarioDimensions
rngState = rng;
rngCleaner = onCleanup(@() rng(rngState));
rng(0);

testData = load('protons_testData.mat');
testData.pln.propDoseCalc.engine = 'MCsquare';
testData.pln.multScen = helper_createAngularRandomScenario(testData.ct, testData.stf);

exception = helper_captureException(@() helper_calcDoseInfluence(testData));
expectedMessagePart = 'does not support active scenario dimension';

assertEqual(exception.identifier, 'matRad:Error');
assertTrue(~isempty(strfind(exception.message, expectedMessagePart)));
assertTrue(~isempty(strfind(exception.message, 'gantry')) || ...
           ~isempty(strfind(exception.message, 'couch')));

function test_nonMigratedEngineRejectsCompactLegacyScenarioStorage
testData = load('protons_testData.mat');
testData.pln.propDoseCalc.engine = 'MCsquare';
testData.pln.multScen = helper_createCompactLegacyScenario(testData.ct);

exception = helper_captureException(@() helper_calcDoseInfluence(testData));
expectedMessagePart = 'supports only legacy ct/setup/range scenario grid storage';

assertEqual(exception.identifier, 'matRad:Error');
assertTrue(~isempty(strfind(exception.message, expectedMessagePart)));

function test_particleImptSteeringRejectsCompactAngularScenarios
rngState = rng;
rngCleaner = onCleanup(@() rng(rngState));
rng(0);

testData = load('protons_testData.mat');
testData.pln.propStf.generator = 'ParticleIMPT';
testData.pln.multScen = helper_createAngularRandomScenario(testData.ct, testData.stf);

exception = helper_captureException(@() matRad_generateStf(testData.ct, testData.cst, testData.pln));

assertEqual(exception.identifier, 'matRad:Error');
assertTrue(~isempty(strfind(exception.message, 'only legacy ct/setup/range')));
assertTrue(~isempty(strfind(exception.message, 'Compact or angular scenario dimensions')));

function dij = helper_calcDoseInfluence(testData)
dij = matRad_calcDoseInfluence(testData.ct, testData.cst, testData.stf, testData.pln);

function multScen = helper_createAngularRandomScenario(ct, stf)
multScen = matRad_RandomScenarios(ct);
multScen.nSamples = 2;
multScen.numOfBeams = numel(stf);
multScen.scenarioDimensionActive = {'ct', 'setup', 'range', 'gantry', 'couch'};

function multScen = helper_createCompactLegacyScenario(ct)
multScen = matRad_NominalScenario(ct);
components = matRad_createScenarioComponents([1 1 1], 1, 1);
scenarioValues = [0 0 0 0 0; 1 0 0 0 0];
ctScenIds = [1; 1];
scenarioWeights = [0.5; 0.5];
scenForProb = [ctScenIds scenarioValues];
storageSubscripts = [(1:2)' ones(2, 1)];
storageSize = [2 1];

multScen.setScenarioRealizations(components, scenarioValues, ctScenIds, ...
                                 scenarioWeights, scenarioWeights, scenForProb, ...
                                 storageSubscripts, storageSize, 'compact-realization');

function exception = helper_captureException(functionHandle)
exception = [];
try
    functionHandle();
catch exception
end

assertFalse(isempty(exception));
