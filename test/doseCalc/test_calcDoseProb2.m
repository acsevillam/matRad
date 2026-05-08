function test_suite = test_calcDoseProb2

test_functions = localfunctions();

initTestSuite;

function test_prob2_physical_dose_expected_and_omega
    [ct,cst,pln,dij,cfg] = singleCtFixture();

    [plnOut,dijProb2Context] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);
    dijProb2 = plnOut.propOpt.dij_prob2;

    expected = [2.5 0; 0 3.5; 1.75 1; 0 4];
    assertElementsAlmostEqual(full(dijProb2.expected),expected,'absolute',1e-12);
    assertElementsAlmostEqual(full(dijProb2.Omega{1}),diag([0.75 0.75]), ...
        'absolute',1e-12);
    assertElementsAlmostEqual(full(dijProb2.Omega{2}), ...
        [0.1875 0; 0 3],'absolute',1e-12);
    assertEqual(dijProb2.quantity,'physicalDose');
    assertEqual(dijProb2.quantityField,'physicalDose');
    assertEqual(dijProb2.probabilisticMode,'PROB2');
    assertElementsAlmostEqual(dijProb2.scenarioWeights,[0.25;0.75], ...
        'absolute',1e-12);
    assertElementsAlmostEqual(full(dijProb2Context.physicalDose{1}), ...
        expected,'absolute',1e-12);
    assertEqual(dijProb2Context.numOfScenarios,1);
    assertTrue(isa(dijProb2Context.scenarioModel,'matRad_NominalScenario'));

function test_prob2_batch_size_does_not_change_result
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.BatchSize = 1;
    [plnBatch,~] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);

    cfg.BatchSize = 99;
    [plnFull,~] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(plnBatch.propOpt.dij_prob2.expected), ...
        full(plnFull.propOpt.dij_prob2.expected),'absolute',1e-12);
    assertElementsAlmostEqual(full(plnBatch.propOpt.dij_prob2.Omega{2}), ...
        full(plnFull.propOpt.dij_prob2.Omega{2}),'absolute',1e-12);

function test_prob2_voi_rows_follow_overlap_priorities
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cst{1,5}.Priority = 1;
    cst{2,5}.Priority = 2;
    cst{1,6} = {struct()};
    cst{2,4}{1} = [2;3;4];
    cst{2,6} = {struct()};

    [plnOut,~] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);
    dijProb2 = plnOut.propOpt.dij_prob2;

    assertEqual(dijProb2.voiSubIx{1},[1;2]);
    assertEqual(dijProb2.voiSubIx{2},[3;4]);

function test_prob2_const_rbe_scales_expected_and_omega
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    dij.RBE = 1.1;
    cfg.Quantity = 'RBExD';

    [plnOut,dijProb2Context] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);
    dijProb2 = plnOut.propOpt.dij_prob2;

    assertElementsAlmostEqual(full(dijProb2.expected(1:2,:)), ...
        1.1*[2.5 0; 0 3.5],'absolute',1e-12);
    assertElementsAlmostEqual(full(dijProb2.Omega{1}), ...
        1.21*diag([0.75 0.75]),'absolute',1e-12);
    assertEqual(dijProb2.quantity,'RBExD');
    assertEqual(dijProb2Context.RBE,1);

function test_prob2_omega_matches_centered_scenario_accumulation
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    dij.totalNumOfBixels = 3;
    dij.physicalDose{1} = sparse([1 0 5; 0 2 7; 1 1 3; 0 1 9]);
    dij.physicalDose{2} = sparse([3 0 5; 0 4 7; 2 1 3; 0 5 9]);

    [plnOut,~] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);
    dijProb2 = plnOut.propOpt.dij_prob2;
    targetRows = cst{1,4}{1};
    expectedRows = dijProb2.expected(targetRows,:);
    scenarioRows = {dij.physicalDose{1}(targetRows,:), ...
        dij.physicalDose{2}(targetRows,:)};
    expectedOmega = manualCenteredOmega(scenarioRows, ...
        dijProb2.scenarioWeights,expectedRows,dij.totalNumOfBixels);

    assertEqual(size(dijProb2.Omega{1}),[3 3]);
    assertElementsAlmostEqual(full(dijProb2.Omega{1}), ...
        full(expectedOmega),'absolute',1e-12);
    assertEqual(nnz(dijProb2.Omega{1}(:,3)),0);
    assertEqual(nnz(dijProb2.Omega{1}(3,:)),0);

function test_prob2_rejects_unsupported_quantity
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.Quantity = 'effect';

    assertExceptionThrown(@() matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function [ct,cst,pln,dij,cfg] = singleCtFixture()
    ct.numOfCtScen = 1;
    ct.cubeDim = [2 2 1];
    ct.resolution = struct('x',1,'y',1,'z',1);
    ct.refScen = 1;

    cst = cell(2,6);
    cst = addStructure(cst,1,'PTV','TARGET',[1;2]);
    cst = addStructure(cst,2,'OAR','OAR',[3;4]);

    pln.bioParam.quantityOpt = 'physicalDose';
    scenarioValues = [0 0 0 0 0; 1 0 0 0 0];
    pln.multScen = fixtureScenarioModel(ct,[1;1],scenarioValues, ...
        [1 1 1; 1 2 1],true(1,2,1),[0.25;0.75]);

    dij = baseDij(ct.cubeDim,2);
    dij.physicalDose = cell(1,2,1);
    dij.physicalDose{1} = sparse([1 0; 0 2; 1 1; 0 1]);
    dij.physicalDose{2} = sparse([3 0; 0 4; 2 1; 0 5]);
    cfg.BatchSize = 2;

function dij = baseDij(dim,numBixels)
    dij.doseGrid.dimensions = dim;
    dij.doseGrid.numOfVoxels = prod(dim);
    dij.doseGrid.resolution = struct('x',1,'y',1,'z',1);
    dij.doseGrid.x = 1:dim(2);
    dij.doseGrid.y = 1:dim(1);
    dij.doseGrid.z = 1:dim(3);
    dij.ctGrid = dij.doseGrid;
    dij.totalNumOfBixels = numBixels;

function cst = addStructure(cst,rowIx,name,type,voxels)
    cst{rowIx,1} = rowIx;
    cst{rowIx,2} = name;
    cst{rowIx,3} = type;
    cst{rowIx,4} = {voxels(:)};
    cst{rowIx,5} = struct();
    cst{rowIx,6} = {};

function scenarioModel = fixtureScenarioModel(ct,ctScenIds,scenarioValues,linearMask,scenMask,scenarioWeights)
    scenarioModel = matRad_NominalScenario(ct);
    dimensions = matRad_createScenarioComponents([1 1 1],1,1);
    scenForProb = [ctScenIds(:) scenarioValues];
    scenarioModel.setScenarioRealizations(dimensions,scenarioValues,ctScenIds, ...
        scenarioWeights(:),scenarioWeights(:),scenForProb,linearMask,scenMask);

function omega = manualCenteredOmega(scenarioRows,scenarioWeights,expectedRows,numBixels)
    omega = sparse(numBixels,numBixels);
    for s = 1:numel(scenarioRows)
        centeredRows = scenarioRows{s} - expectedRows;
        omega = omega + scenarioWeights(s).*(centeredRows' * centeredRows);
    end
    omega = sparse(0.5.*(omega + omega'));
