function test_suite = test_scenarioModel

test_functions=localfunctions();

initTestSuite;

%Add automated instance tests for avaiable models
ct.numOfCtScen = 5;
models = {matRad_NominalScenario(),matRad_WorstCaseScenarios(),matRad_ImportanceScenarios(), ...
    matRad_TruncatedImportanceScenarios(),matRad_RandomScenarios(),matRad_NominalScenario(ct), ...
    matRad_WorstCaseScenarios(ct),matRad_ImportanceScenarios(ct), ...
    matRad_TruncatedImportanceScenarios(ct),matRad_RandomScenarios(ct)};

instanceTests = cellfun(@func2str,test_functions,'UniformOutput',false);
funIx = ~cellfun(@isempty,strfind(instanceTests,'instanceTest_'));
instanceTests = instanceTests(funIx);
instanceTestHandles = cellfun(@str2func,instanceTests,'UniformOutput',false);

for m = 1:numel(models)
    model = models{m};
    for i = 1:numel(instanceTests)
        funHandle = @() instanceTestHandles{i}(model);
        test_case=MOxUnitFunctionHandleTestCase([instanceTests{i} '_' class(model)],mfilename,funHandle);
        test_suite=addTest(test_suite, test_case);
    end
end

function assignmentTestHelper(model,property,value)
    model.(property) = value;

function stf = scenarioStfSignature(numOfBeams)
    if nargin < 1
        numOfBeams = 1;
    end

    stf = struct('isoCenter',{},'gantryAngle',{},'couchAngle',{});
    for i = 1:numOfBeams
        stf(i).isoCenter = [0 0 0];
        stf(i).gantryAngle = 0;
        stf(i).couchAngle = 0;
    end

function key = scenarioStfKey(stf)
    values = [];
    for i = 1:numel(stf)
        values = [values, stf(i).isoCenter, stf(i).gantryAngle, stf(i).couchAngle];
    end
    key = sprintf('%.10g,',values);

function test_scenarioAbstract
    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@() matRad_ScenarioModel(),'');
    else
        assertExceptionThrown(@() matRad_ScenarioModel(),'MATLAB:class:abstract');
    end

function test_scenarioAbstractAvailableTypes()
    availableTypes = {'nomScen','wcScen','impScen','truncatedImpScen','rndScen'};
    assertTrue(iscell(availableTypes))
    assertTrue(all(cellfun(@ischar,availableTypes)));

    for i = 1:numel(availableTypes)
        model = matRad_createScenarioModel([],availableTypes{i});
        assertTrue(isa(model,'matRad_ScenarioModel'));
        assertEqual(model.name,availableTypes{i});
    end

function test_supportedScenarioDimensions_includeAngularButDefaultsDoNot
    assertEqual(matRad_defaultScenarioUncertaintyDimensions(), ...
        {'ct','setup','range','gantry','couch'});
    assertEqual(matRad_defaultScenarioActiveDimensions(), ...
        {'ct','setup','range'});

    model = matRad_RandomScenarios();
    assertEqual(model.scenarioDimensionActive,{'ct','setup','range'});

function test_griddedScenariosRejectAngularDimensions
    model = matRad_WorstCaseScenarios();
    model.numOfBeams = 2;

    assertExceptionThrown(@() assignmentTestHelper(model,'scenarioDimensionActive', ...
        {'ct','setup','range','gantry'}),'matRad:Error');

function test_extractSingleScenario_accepts_sparse_ct_scenario_probabilities
    ct.numOfCtScen = 3;
    model = matRad_RandomScenarios(ct);
    model.ctScenProb = [2 1];

    scenario = model.extractSingleScenario(1);

    assertEqual(scenario.ctScenProb,[2 1]);
    assertEqual(scenario.scenarioCtScenIds,2);
    assertEqual(scenario.sub2scenIx(1,1,1),2);
    assertEqual(scenario.sub2scenIx(1,1,1,'position'),2);
    assertEqual(scenario.sub2scenIx(2,1,1,'id'),2);
    assertExceptionThrown(@() scenario.sub2scenIx(2,1,1,'position'),'matRad:Error');

function test_omittedCtScenariosAreInactive
    ct.numOfCtScen = 3;
    availableTypes = {'nomScen','wcScen','impScen','truncatedImpScen','rndScen'};

    for i = 1:numel(availableTypes)
        model = matRad_createScenarioModel(ct,availableTypes{i});
        model.ctScenProb = [1 0.25; 3 0.75];

        inactiveCtMask = model.scenMask(2,:,:);
        assertFalse(any(model.scenarioCtScenIds == 2));
        assertFalse(any(model.linearMask(:,1) == 2));
        assertFalse(any(inactiveCtMask(:)));
        assertEqual(model.totNumScen,size(model.scenarioValues,1));
        assertEqual(model.totNumScen,numel(model.scenWeight));
        assertTrue(all(model.scenWeight > 0));
        assertElementsAlmostEqual(sum(model.scenWeight),1);
    end

function test_inactiveScenarioDimensionsAllowZeroScale
    model = matRad_ImportanceScenarios();
    model.scenarioDimensionActive = {'ct','setup'};
    model.rangeAbsSD = 0;
    model.rangeRelSD = 0;

    names = {model.scenarioComponents.name};
    rangeAbsIx = find(strcmp(names,'range.absolute'),1,'first');
    rangeRelIx = find(strcmp(names,'range.relative'),1,'first');

    assertFalse(model.scenarioComponents(rangeAbsIx).active);
    assertFalse(model.scenarioComponents(rangeRelIx).active);
    assertEqual(model.totNumRangeScen,1);
    assertTrue(all(model.absRangeShift == 0));
    assertTrue(all(model.relRangeShift == 0));
    assertElementsAlmostEqual(sum(model.scenWeight),1);

    model.scenarioDimensionActive = {};
    model.shiftSD = [0 0 0];
    assertFalse(any([model.scenarioComponents.active]));
    assertTrue(all(model.isoShift(:) == 0));
    assertTrue(all(model.absRangeShift == 0));
    assertTrue(all(model.relRangeShift == 0));
    assertElementsAlmostEqual(sum(model.scenWeight),1);

function test_activeScenarioDimensionsRejectZeroScale
    model = matRad_ImportanceScenarios();

    assertExceptionThrown(@() assignmentTestHelper(model,'rangeAbsSD',0),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeRelSD',0),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'shiftSD',[0 2.25 2.25]),'matRad:Error');

function test_activatingZeroScaleDimensionFails
    model = matRad_ImportanceScenarios();
    model.scenarioDimensionActive = {'ct','setup'};
    model.rangeRelSD = 0;

    assertExceptionThrown(@() assignmentTestHelper(model,'scenarioDimensionActive', ...
        {'ct','setup','range'}),'matRad:Error');

function test_singleCtSetupScenarioSubscriptsUseSetupDimension
    model = matRad_WorstCaseScenarios();
    model.scenarioDimensionActive = {'ct','setup'};
    model.shiftSD = [1 2 3];
    model.wcSigma = 1;

    assertEqual(size(model.scenMask),[1 7]);
    assertEqual(model.totNumShiftScen,7);
    assertEqual(model.totNumRangeScen,1);

    fullScenIx = zeros(model.totNumShiftScen,1);
    for setupScenIx = 1:model.totNumShiftScen
        fullScenIx(setupScenIx) = model.getDijScenarioIndexBySubscripts( ...
            1,setupScenIx,1,'position');
    end

    assertEqual(fullScenIx,(1:model.totNumShiftScen)');

function test_applyScenarioToStfIgnoresCtAndRange
    ct.numOfCtScen = 2;
    model = matRad_RandomScenarios(ct);
    model.nSamples = 3;
    model.randomSeed = 17;
    model.scenarioDimensionActive = {'ct','range'};

    stf = scenarioStfSignature();
    scenarioIds = model.scenarioIds();
    stfKeys = cell(numel(scenarioIds),1);
    for i = 1:numel(scenarioIds)
        stfKeys{i} = scenarioStfKey(model.applyScenarioToStf(scenarioIds(i),stf));
    end

    assertEqual(numel(unique(stfKeys)),1);

function test_applyScenarioToStfFollowsSetupNotCtRangeProduct
    ct.numOfCtScen = 2;
    model = matRad_RandomScenarios(ct);
    model.nSamples = 4;
    model.randomSeed = 19;
    model.scenarioDimensionActive = {'ct','setup','range'};

    stf = scenarioStfSignature();
    scenarioIds = model.scenarioIds();
    stfKeys = cell(numel(scenarioIds),1);
    for i = 1:numel(scenarioIds)
        stfKeys{i} = scenarioStfKey(model.applyScenarioToStf(scenarioIds(i),stf));
    end

    uniqueKeys = unique(stfKeys);
    assertEqual(numel(uniqueKeys),model.nSamples);
    for i = 1:numel(uniqueKeys)
        assertEqual(sum(strcmp(stfKeys,uniqueKeys{i})),ct.numOfCtScen);
    end

function test_nominalAngularScenarioUsesNominalStfSignature
    nominalModel = matRad_NominalScenario();
    nominalModel.numOfBeams = 2;

    angularModel = matRad_RandomScenarios();
    angularModel.numOfBeams = 2;
    angularModel.nSamples = 3;
    angularModel.randomSeed = 23;
    angularModel.scenarioDimensionActive = {'ct','gantry','couch'};
    angularModel.includeNominalScenario = true;

    nominalScenarioIds = angularModel.getNominalScenarioIds();

    assertFalse(isempty(nominalScenarioIds));
    stf = scenarioStfSignature(2);
    assertEqual(scenarioStfKey(angularModel.applyScenarioToStf(nominalScenarioIds(1),stf)), ...
        scenarioStfKey(nominalModel.applyScenarioToStf(1,stf)));

function instanceTest_listAllScenarios(model)
    model.listAllScenarios();
    assertTrue(true); %Will be reached if above call does not fail

function instanceTest_relRangeUncertainty(model)
    newValue = 0.01;
    model.rangeRelSD = newValue;
    getValue = model.rangeRelSD;
    assertEqual(newValue,getValue);
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeRelSD','a'),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeRelSD',ones(2)),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeRelSD',-1),'matRad:Error');

function instanceTest_absRangeUncertainty(model)
    newValue = 2;
    model.rangeAbsSD = newValue;
    getValue = model.rangeAbsSD;
    assertEqual(newValue,getValue);
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeAbsSD','a'),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeAbsSD',ones(2)),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeAbsSD',-1),'matRad:Error');

function instanceTest_shiftUncertainty(model)
    newValue = [1 1 1];
    model.shiftSD = newValue;
    getValue = model.shiftSD;
    assertEqual(newValue,getValue);
    assertExceptionThrown(@() assignmentTestHelper(model,'shiftSD','a'),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'shiftSD',ones(3,3)),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'shiftSD',ones(3,1)),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'shiftSD',[-1 2 2]),'matRad:Error');

function instanceTest_wcSigma(model)
    newValue = 5;
    model.wcSigma = newValue;
    getValue = model.wcSigma;
    assertEqual(newValue,getValue);
    assertExceptionThrown(@() assignmentTestHelper(model,'wcSigma','a'),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'wcSigma',-1),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'wcSigma',ones(2)),'matRad:Error');

function instanceTest_ctScenProb(model)
    newValue = (1:model.numOfCtScen)';
    newValue = newValue ./ sum(newValue);
    model.ctScenProb(:,2) = newValue;
    getValue = model.ctScenProb(:,2);
    assertEqual(newValue,getValue);
    assertExceptionThrown(@() assignmentTestHelper(model,'ctScenProb','a'),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'ctScenProb',ones(1,5)),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'ctScenProb',-1*ones(model.numOfCtScen,1)),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'ctScenProb',[(1:model.numOfCtScen)' zeros(model.numOfCtScen,1)]),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'ctScenProb',[1 0.5; 2 0]),'matRad:Error');

function instanceTest_TYPE(model)
    %assertWarning(@() model.TYPE,'matRad:Deprecated');
    assertEqual(model.TYPE,model.name);

function instanceTest_wcFactor(model)
    %assertWarning(@() model.wcFactor,'matRad:Deprecated');
    assertEqual(model.TYPE,model.name);
