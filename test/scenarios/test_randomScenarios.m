function test_suite = test_randomScenarios

test_functions=localfunctions();

initTestSuite;

function test_randomScenarioConstructor
    scenario = matRad_RandomScenarios();

    %Defaults
    nSamples = 10;

    assertTrue(isa(scenario, 'matRad_RandomScenarios'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'rndScen');
    %Test correct standard values & sizes
    assertEqual(scenario.ctScenProb, [1 1]);
    assertEqual(scenario.numOfCtScen, 1);
    assertEqual(scenario.totNumScen, nSamples);
    assertEqual(scenario.totNumShiftScen, nSamples);
    assertEqual(scenario.totNumRangeScen, nSamples);
    assertEqual(size(scenario.relRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.absRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.isoShift),[scenario.totNumScen,3]);
    %assertEqual(scenario.isoShift, 2.25 * ones(1,3));
    assertEqual(scenario.maxAbsRangeShift, max(scenario.absRangeShift));
    assertEqual(scenario.maxRelRangeShift, max(scenario.relRangeShift));
    assertEqual(size(scenario.scenMask), [scenario.numOfCtScen,scenario.totNumShiftScen,scenario.totNumRangeScen]);
    %assertEqual(scenario.scenMask, true(1,1,1));
    assertEqual(size(scenario.linearMask), [scenario.totNumScen,3]);
    assertScenarioMaskDerivedFromLinearMask(scenario);

    %assertEqual(ind2sub(find()))
    %assertEqual(scenario.linearMask, [1 1 1]);
    assertElementsAlmostEqual(scenario.scenProb,helper_mvarGauss(scenario));

    tmp = [scenario.scenarioCtScenIds scenario.isoShift scenario.absRangeShift scenario.relRangeShift];
    assertEqual(scenario.scenForProb,tmp);

    assertEqual(numel(unique(scenario.scenWeight)),scenario.numOfCtScen);

function test_randomScenarioConstructorWithCt
    n = 5;
    ct = struct('numOfCtScen',n);

    %Defaults
    nSamples = 10;

    scenario = matRad_RandomScenarios(ct);
    assertTrue(isa(scenario, 'matRad_RandomScenarios'));
    assertTrue(isa(scenario, 'matRad_ScenarioModel'));
    assertEqual(scenario.name, 'rndScen');
    %Test correct standard values & sizes
    assertEqual(scenario.ctScenProb, [(1:n)' ones(n,1)./n]);
    assertEqual(scenario.numOfCtScen, n);
    assertEqual(scenario.totNumScen, nSamples*n);
    assertEqual(scenario.totNumShiftScen, nSamples);
    assertEqual(scenario.totNumRangeScen, nSamples);
    assertEqual(size(scenario.relRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.absRangeShift),[scenario.totNumScen,1]);
    assertEqual(size(scenario.isoShift),[scenario.totNumScen,3]);
    %assertEqual(scenario.isoShift, 2.25 * ones(1,3));
    assertEqual(scenario.maxAbsRangeShift, max(scenario.absRangeShift));
    assertEqual(scenario.maxRelRangeShift, max(scenario.relRangeShift));
    assertEqual(size(scenario.scenMask), [scenario.numOfCtScen,scenario.totNumShiftScen,scenario.totNumRangeScen]);
    %assertEqual(scenario.scenMask, true(1,1,1));
    assertEqual(size(scenario.linearMask), [scenario.totNumScen,3]);
    assertScenarioMaskDerivedFromLinearMask(scenario);

    %assertEqual(ind2sub(find()))
    %assertEqual(scenario.linearMask, [1 1 1]);
    assertElementsAlmostEqual(scenario.scenProb,helper_mvarGauss(scenario));

    tmp = [scenario.scenarioCtScenIds scenario.isoShift scenario.absRangeShift scenario.relRangeShift];
    assertEqual(scenario.scenForProb,tmp);
    assertEqual(numel(unique(scenario.scenWeight)),numel(unique(scenario.ctScenProb(:,2))));

function test_randomScenarioSetupOnlySingleCtLinearMask
    nSamples = 4;
    scenario = matRad_RandomScenarios();
    scenario.nSamples = nSamples;
    scenario.randomSeed = 31;
    scenario.scenarioDimensionActive = {'ct','setup'};

    expectedLinearMask = [ones(nSamples,1) (1:nSamples)' ones(nSamples,1)];

    assertEqual(scenario.totNumShiftScen,nSamples);
    assertEqual(scenario.totNumRangeScen,1);
    assertEqual(scenario.linearMask,expectedLinearMask);
    assertEqual(scenario.scenarioCtScenIds,ones(nSamples,1));
    assertScenarioMaskDerivedFromLinearMask(scenario);

function test_randomScenarioRangeOnlySingleCtLinearMask
    nSamples = 4;
    scenario = matRad_RandomScenarios();
    scenario.nSamples = nSamples;
    scenario.randomSeed = 37;
    scenario.scenarioDimensionActive = {'ct','range'};

    expectedLinearMask = [ones(nSamples,1) ones(nSamples,1) (1:nSamples)'];

    assertEqual(scenario.totNumShiftScen,1);
    assertEqual(scenario.totNumRangeScen,nSamples);
    assertEqual(scenario.linearMask,expectedLinearMask);
    assertEqual(scenario.scenarioCtScenIds,ones(nSamples,1));
    assertScenarioMaskDerivedFromLinearMask(scenario);

function test_randomScenarioSetupRangeSingleCtLinearMask
    nSamples = 4;
    scenario = matRad_RandomScenarios();
    scenario.nSamples = nSamples;
    scenario.randomSeed = 41;
    scenario.scenarioDimensionActive = {'ct','setup','range'};

    expectedLinearMask = [ones(nSamples,1) (1:nSamples)' (1:nSamples)'];

    assertEqual(scenario.totNumShiftScen,nSamples);
    assertEqual(scenario.totNumRangeScen,nSamples);
    assertEqual(scenario.linearMask,expectedLinearMask);
    assertEqual(scenario.scenarioCtScenIds,ones(nSamples,1));
    assertScenarioMaskDerivedFromLinearMask(scenario);

function test_randomScenarioSetupOnlyMultiCtLinearMask
    nCtScen = 3;
    nSamples = 3;
    ct = struct('numOfCtScen',nCtScen);
    scenario = matRad_RandomScenarios(ct);
    scenario.nSamples = nSamples;
    scenario.randomSeed = 43;
    scenario.scenarioDimensionActive = {'ct','setup'};

    expectedCtScenIds = kron((1:nCtScen)',ones(nSamples,1));
    expectedLinearMask = [expectedCtScenIds repmat((1:nSamples)',nCtScen,1) ...
        ones(nSamples*nCtScen,1)];

    assertEqual(scenario.totNumShiftScen,nSamples);
    assertEqual(scenario.totNumRangeScen,1);
    assertEqual(scenario.linearMask,expectedLinearMask);
    assertEqual(scenario.scenarioCtScenIds,expectedCtScenIds);
    assertScenarioMaskDerivedFromLinearMask(scenario);

function test_randomScenarioSparseCtSetupOnlyLinearMask
    nSamples = 3;
    ct = struct('numOfCtScen',3);
    scenario = matRad_RandomScenarios(ct);
    scenario.nSamples = nSamples;
    scenario.randomSeed = 47;
    scenario.ctScenProb = [2 1];
    scenario.scenarioDimensionActive = {'ct','setup'};

    expectedLinearMask = [2*ones(nSamples,1) (1:nSamples)' ones(nSamples,1)];

    assertEqual(scenario.numOfCtScen,1);
    assertEqual(scenario.linearMask,expectedLinearMask);
    assertEqual(scenario.scenarioCtScenIds,2*ones(nSamples,1));
    assertFalse(any(scenario.linearMask(:,1) == 1));
    assertFalse(any(scenario.linearMask(:,1) == 3));
    assertScenarioMaskDerivedFromLinearMask(scenario);

function test_randomScenarioCtOnlyUsesOneScenarioPerCt
    nCtScen = 3;
    ct = struct('numOfCtScen',nCtScen);
    scenario = matRad_RandomScenarios(ct);
    scenario.nSamples = 4;
    scenario.scenarioDimensionActive = {'ct'};

    expectedLinearMask = [(1:nCtScen)' ones(nCtScen,2)];

    assertEqual(scenario.totNumShiftScen,1);
    assertEqual(scenario.totNumRangeScen,1);
    assertEqual(scenario.totNumScen,nCtScen);
    assertEqual(scenario.linearMask,expectedLinearMask);
    assertEqual(scenario.scenarioCtScenIds,(1:nCtScen)');
    assertElementsAlmostEqual(scenario.scenWeight,ones(nCtScen,1)./nCtScen);
    assertScenarioMaskDerivedFromLinearMask(scenario);

function test_randomScenarioAngularComponents
    nSamples = 4;
    scenario = matRad_RandomScenarios();
    scenario.numOfBeams = 2;
    scenario.nSamples = nSamples;
    scenario.randomSeed = 42;
    scenario.scenarioDimensionActive = {'ct','gantry','couch'};

    names = {scenario.scenarioComponents.name};
    assertTrue(any(strcmp(names,'gantry.beam1')));
    assertTrue(any(strcmp(names,'gantry.beam2')));
    assertTrue(any(strcmp(names,'couch.beam1')));
    assertTrue(any(strcmp(names,'couch.beam2')));
    assertEqual(scenario.totNumShiftScen,1);
    assertEqual(scenario.totNumRangeScen,1);
    assertEqual(scenario.totNumGantryScen,nSamples);
    assertEqual(scenario.totNumCouchScen,nSamples);
    assertEqual(scenario.totNumScen,nSamples);
    assertEqual(size(scenario.linearMask),[nSamples 5]);
    assertEqual(scenario.linearMask(:,2:3),ones(nSamples,2));
    assertEqual(scenario.linearMask(:,4),(1:nSamples)');
    assertEqual(scenario.linearMask(:,5),(1:nSamples)');
    assertEqual(size(scenario.gantryAngleOffset),[nSamples 2]);
    assertEqual(size(scenario.couchAngleOffset),[nSamples 2]);
    assertScenarioMaskDerivedFromLinearMask(scenario);

function test_randomScenarioGantryOnlyAngularComponents
    nSamples = 4;
    scenario = matRad_RandomScenarios();
    scenario.numOfBeams = 2;
    scenario.nSamples = nSamples;
    scenario.randomSeed = 53;
    scenario.scenarioDimensionActive = {'ct','gantry'};

    assertEqual(scenario.totNumShiftScen,1);
    assertEqual(scenario.totNumRangeScen,1);
    assertEqual(scenario.totNumGantryScen,nSamples);
    assertEqual(scenario.totNumCouchScen,1);
    assertEqual(scenario.totNumScen,nSamples);
    assertEqual(size(scenario.linearMask),[nSamples 5]);
    assertEqual(scenario.linearMask(:,2:3),ones(nSamples,2));
    assertEqual(scenario.linearMask(:,4),(1:nSamples)');
    assertEqual(scenario.linearMask(:,5),ones(nSamples,1));
    assertEqual(size(scenario.gantryAngleOffset),[nSamples 2]);
    assertEqual(scenario.couchAngleOffset,zeros(nSamples,2));
    assertScenarioMaskDerivedFromLinearMask(scenario);

function test_randomScenarioCouchOnlyAngularComponents
    nSamples = 4;
    scenario = matRad_RandomScenarios();
    scenario.numOfBeams = 2;
    scenario.nSamples = nSamples;
    scenario.randomSeed = 59;
    scenario.scenarioDimensionActive = {'ct','couch'};

    assertEqual(scenario.totNumShiftScen,1);
    assertEqual(scenario.totNumRangeScen,1);
    assertEqual(scenario.totNumGantryScen,1);
    assertEqual(scenario.totNumCouchScen,nSamples);
    assertEqual(scenario.totNumScen,nSamples);
    assertEqual(size(scenario.linearMask),[nSamples 5]);
    assertEqual(scenario.linearMask(:,2:4),ones(nSamples,3));
    assertEqual(scenario.linearMask(:,5),(1:nSamples)');
    assertEqual(scenario.gantryAngleOffset,zeros(nSamples,2));
    assertEqual(size(scenario.couchAngleOffset),[nSamples 2]);
    assertScenarioMaskDerivedFromLinearMask(scenario);

function test_randomScenarioAngularActivationCanWaitForBeamCount
    scenario = matRad_RandomScenarios();
    scenario.scenarioDimensionActive = {'ct','gantry','couch'};
    scenario.numOfBeams = 2;

    assertEqual(scenario.scenarioDimensionActive,{'ct','gantry','couch'});
    assertTrue(any(strcmp({scenario.scenarioComponents.name},'gantry.beam1')));
    assertTrue(any(strcmp({scenario.scenarioComponents.name},'couch.beam2')));

function test_randomScenarioAngularSeedReproducibility
    scenarioA = matRad_RandomScenarios();
    scenarioA.numOfBeams = 2;
    scenarioA.nSamples = 5;
    scenarioA.randomSeed = 7;
    scenarioA.scenarioDimensionActive = {'ct','gantry','couch'};

    scenarioB = matRad_RandomScenarios();
    scenarioB.numOfBeams = 2;
    scenarioB.nSamples = 5;
    scenarioB.randomSeed = 7;
    scenarioB.scenarioDimensionActive = {'ct','gantry','couch'};

    assertElementsAlmostEqual(scenarioA.scenarioValues,scenarioB.scenarioValues);

function test_randomScenarioAngularIncludeNominalScenario
    scenario = matRad_RandomScenarios();
    scenario.numOfBeams = 2;
    scenario.nSamples = 5;
    scenario.randomSeed = 11;
    scenario.scenarioDimensionActive = {'ct','gantry','couch'};
    scenario.includeNominalScenario = true;

    assertEqual(scenario.scenarioValues(1,:),zeros(1,size(scenario.scenarioValues,2)));
    assertEqual(scenario.gantryAngleOffset(1,:),[0 0]);
    assertEqual(scenario.couchAngleOffset(1,:),[0 0]);

function test_randomScenarioExtractSingleAngularScenario
    scenario = matRad_RandomScenarios();
    scenario.numOfBeams = 2;
    scenario.nSamples = 4;
    scenario.randomSeed = 21;
    scenario.scenarioDimensionActive = {'ct','gantry','couch'};

    scenarioId = 2;
    singleScenario = scenario.extractSingleScenario(scenarioId);
    ctScenId = scenario.scenarioCtScenIds(scenarioId);
    stf = randomScenarioTestStf(scenario.numOfBeams);
    originalStf = scenario.applyScenarioToStf(scenarioId,stf);
    extractedStf = singleScenario.applyScenarioToStf(1,stf);

    assertEqual(singleScenario.scenarioDimensionActive,scenario.scenarioDimensionActive);
    assertEqual(singleScenario.numOfBeams,scenario.numOfBeams);
    assertEqual(singleScenario.gantryAngleOffset,scenario.gantryAngleOffset(scenarioId,:));
    assertEqual(singleScenario.couchAngleOffset,scenario.couchAngleOffset(scenarioId,:));
    assertEqual(singleScenario.linearMask,[ctScenId 1 1 1 1]);
    assertScenarioMaskDerivedFromLinearMask(singleScenario);
    assertElementsAlmostEqual(randomScenarioStfSignature(extractedStf), ...
        randomScenarioStfSignature(originalStf));


function test_randomScenarioExtractSingleScenario
    refScen = matRad_RandomScenarios();
    for scenNum = 1:refScen.totNumScen
        scenario = refScen.extractSingleScenario(scenNum);
        assertTrue(isa(scenario, 'matRad_NominalScenario'));
        ctScenId = refScen.scenarioCtScenIds(scenNum);
        ctScenProbIx = find(ctScenId == refScen.ctScenProb(:,1));
        assertEqual(scenario.ctScenProb, refScen.ctScenProb(ctScenProbIx,:));
        assertEqual(scenario.numOfCtScen, 1);
        assertEqual(scenario.totNumScen, 1);
        assertEqual(scenario.totNumShiftScen, 1);
        assertEqual(scenario.totNumRangeScen, 1);
        assertEqual(scenario.relRangeShift, refScen.relRangeShift(scenNum));
        assertEqual(scenario.absRangeShift, refScen.absRangeShift(scenNum));
        assertEqual(scenario.isoShift, refScen.isoShift(scenNum,:));
        assertEqual(scenario.maxAbsRangeShift, max(abs(refScen.absRangeShift(scenNum))));
        assertEqual(scenario.maxRelRangeShift, max(abs(refScen.relRangeShift(scenNum))));
        assertTrue(scenario.scenMask(ctScenId,1,1));
        assertTrue(numel(find(scenario.scenMask)) == 1);
        assertEqual(scenario.linearMask, [ctScenId 1 1]);
        assertElementsAlmostEqual(scenario.scenProb,helper_mvarGauss(scenario));
        assertEqual(scenario.scenForProb,refScen.scenForProb(scenNum,:));
        assertEqual(scenario.scenWeight, refScen.scenWeight(scenNum));
    end

function test_randomScenarioExtractSingleScenarioWithCtScen
    n = 5;
    ct = struct('numOfCtScen',n);
    refScen = matRad_RandomScenarios(ct);
    for scenNum = 1:refScen.totNumScen
        scenario = refScen.extractSingleScenario(scenNum);
        assertTrue(isa(scenario, 'matRad_NominalScenario'));
        ctScenId = refScen.scenarioCtScenIds(scenNum);
        ctScenProbIx = find(ctScenId == refScen.ctScenProb(:,1));
        assertEqual(scenario.ctScenProb, refScen.ctScenProb(ctScenProbIx,:));
        assertEqual(scenario.numOfCtScen, 1);
        assertEqual(scenario.totNumScen, 1);
        assertEqual(scenario.totNumShiftScen, 1);
        assertEqual(scenario.totNumRangeScen, 1);
        assertEqual(scenario.relRangeShift, refScen.relRangeShift(scenNum));
        assertEqual(scenario.absRangeShift, refScen.absRangeShift(scenNum));
        assertEqual(scenario.isoShift, refScen.isoShift(scenNum,:));
        assertEqual(scenario.maxAbsRangeShift, max(abs(refScen.absRangeShift(scenNum))));
        assertEqual(scenario.maxRelRangeShift, max(abs(refScen.relRangeShift(scenNum))));
        assertTrue(scenario.scenMask(ctScenId,1,1));
        assertTrue(numel(find(scenario.scenMask)) == 1);
        assertEqual(scenario.linearMask, [ctScenId 1 1]);
        assertElementsAlmostEqual(scenario.scenProb,helper_mvarGauss(scenario));
        assertEqual(scenario.scenForProb,refScen.scenForProb(scenNum,:));
        assertEqual(scenario.scenWeight, refScen.scenWeight(scenNum));
    end

function assertScenarioMaskDerivedFromLinearMask(scenario)
    maskSize = paddedScenarioMaskSize(scenario.scenMask,size(scenario.linearMask,2));
    expectedMask = false(maskSize);
    maskSubscripts = mat2cell(scenario.linearMask,size(scenario.linearMask,1), ...
        ones(1,size(scenario.linearMask,2)));
    expectedMask(sub2ind(maskSize,maskSubscripts{:})) = true;

    assertEqual(expectedMask(:),scenario.scenMask(:));
    assertEqual(numel(find(scenario.scenMask)),size(scenario.linearMask,1));

function maskSize = paddedScenarioMaskSize(scenMask,numDims)
    maskSize = size(scenMask);
    if numel(maskSize) < numDims
        maskSize(end+1:numDims) = 1;
    end

function stf = randomScenarioTestStf(numOfBeams)
    stf = struct('isoCenter',{},'gantryAngle',{},'couchAngle',{});
    for b = 1:numOfBeams
        stf(b).isoCenter = [0 0 0];
        stf(b).gantryAngle = 10*b;
        stf(b).couchAngle = 2*b;
    end

function signature = randomScenarioStfSignature(stf)
    signature = [];
    for b = 1:numel(stf)
        signature = [signature stf(b).isoCenter stf(b).gantryAngle stf(b).couchAngle];
    end
