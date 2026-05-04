function test_suite = test_angularRandomScenarios

test_functions=localfunctions();

initTestSuite;

function test_photonPencilBeamDijWithAngularRandomScenarios
    [ct,cst,pln,stf] = loadTestData('photons_testData.mat');
    pln.multScen = angularRandomScenario(ct);
    pln.propDoseCalc.engine = 'SVDPB';

    dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

    assertAngularDij(dij,pln.multScen);

function test_protonPencilBeamDijWithAngularRandomScenarios
    [ct,cst,pln,~] = loadTestData('protons_testData.mat');
    pln.multScen = angularRandomScenario(ct);
    pln.propDoseCalc.engine = 'HongPB';
    stf = matRad_generateStf(ct,cst,pln);

    dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

    assertAngularDij(dij,pln.multScen);

function test_photonPencilBeamDijWithRangeRandom
    [ct,cst,pln,stf] = loadTestData('photons_testData.mat');
    pln.multScen = rangeRandomScenario(ct);
    pln.propDoseCalc.engine = 'SVDPB';

    dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

    assertEqual(numel(dij.physicalDose),numel(dij.scenarioModel.scenMask));
    assertEqual(numel(dij.beamNum),dij.totalNumOfBixels);
    assertTrue(all(isfinite(dij.beamNum)));

function test_monteCarloEnginesRejectAngularRandomScenarios
    [ct,cst,pln,stf] = loadTestData('protons_testData.mat');
    pln.multScen = angularRandomScenario(ct);
    pln.propDoseCalc.engine = 'MCsquare';

    assertExceptionThrown(@() matRad_calcDoseInfluence(ct,cst,stf,pln),'matRad:Error');

function [ct,cst,pln,stf] = loadTestData(fileName)
    testDataPath = fullfile(fileparts(mfilename('fullpath')),'..','testData',fileName);
    data = load(testDataPath,'ct','cst','pln','stf');
    ct = data.ct;
    cst = data.cst;
    pln = data.pln;
    stf = data.stf;

function scenario = angularRandomScenario(ct)
    scenario = matRad_RandomScenarios(ct);
    scenario.nSamples = 2;
    scenario.randomSeed = 13;
    scenario.gantryAngleSD = 0.25;
    scenario.couchAngleSD = 0.25;
    scenario.scenarioDimensionActive = {'ct','gantry','couch'};
    assertEqual(scenario.numOfBeams,0);

function scenario = rangeRandomScenario(ct)
    scenario = matRad_RandomScenarios(ct);
    scenario.nSamples = 2;
    scenario.randomSeed = 29;
    scenario.scenarioDimensionActive = {'ct','range'};

function assertAngularDij(dij,scenario)
    assertEqual(dij.numOfScenarios,scenario.totNumScen);
    assertEqual(dij.scenarioModel.scenarioDimensionActive,scenario.scenarioDimensionActive);
    assertEqual(size(dij.scenarioModel.scenarioValues),size(scenario.scenarioValues));
    assertEqual(size(dij.scenarioModel.linearMask,2),5);
    assertEqual(numel(find(dij.scenarioModel.scenMask)),scenario.totNumScen);
    assertEqual(numel(dij.physicalDose),numel(dij.scenarioModel.scenMask));
    assertEqual(numel(dij.beamNum),dij.totalNumOfBixels);
    assertTrue(all(isfinite(dij.beamNum)));

    scenarioIds = dij.scenarioModel.scenarioIds();
    for i = 1:numel(scenarioIds)
        fullScenIx = dij.scenarioModel.getDijScenarioIndex(scenarioIds(i));
        assertFalse(isempty(dij.physicalDose{fullScenIx}));
        assertEqual(size(dij.physicalDose{fullScenIx},2),dij.totalNumOfBixels);
    end
