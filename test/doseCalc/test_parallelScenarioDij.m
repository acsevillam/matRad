function test_suite = test_parallelScenarioDij

test_functions = localfunctions();

initTestSuite;

function test_assembler_inserts_matrices_at_original_dij_indices
    scenarioModel = fixtureScenarioModel();
    dij1 = scenarioDij([1 0; 0 1]);
    dij2 = scenarioDij([2 0; 0 2]);

    dij = matRad_assembleParallelScenarioDij({dij1;dij2},[1;2],scenarioModel);

    assertEqual(size(dij.physicalDose),size(scenarioModel.scenMask));
    assertElementsAlmostEqual(full(dij.physicalDose{1}), ...
        [1 0; 0 1],'absolute',1e-12);
    assertTrue(isempty(dij.physicalDose{2}));
    assertElementsAlmostEqual(full(dij.physicalDose{3}), ...
        [2 0; 0 2],'absolute',1e-12);
    assertEqual(dij.numOfScenarios,2);
    assertEqual(dij.scenarioIds,[1;2]);
    assertEqual(dij.scenarioModel.scenMask,scenarioModel.scenMask);
    assertEqual(dij.beamNum,dij1.beamNum);

function test_assembler_handles_additional_dose_matrix_fields
    scenarioModel = fixtureScenarioModel();
    dij1 = scenarioDij([1 0; 0 1]);
    dij2 = scenarioDij([2 0; 0 2]);
    dij1.mLETDose = {6*dij1.physicalDose{1}};
    dij2.mLETDose = {7*dij2.physicalDose{1}};
    dij1.mAlphaDose = {2*dij1.physicalDose{1}};
    dij2.mAlphaDose = {3*dij2.physicalDose{1}};
    dij1.mSqrtBetaDose = {4*dij1.physicalDose{1}};
    dij2.mSqrtBetaDose = {5*dij2.physicalDose{1}};
    dij1.customQuantity = {8*dij1.physicalDose{1}};
    dij2.customQuantity = {9*dij2.physicalDose{1}};

    dij = matRad_assembleParallelScenarioDij({dij1;dij2},[1;2],scenarioModel);

    assertElementsAlmostEqual(full(dij.mLETDose{1}), ...
        6*[1 0; 0 1],'absolute',1e-12);
    assertElementsAlmostEqual(full(dij.mLETDose{3}), ...
        7*[2 0; 0 2],'absolute',1e-12);
    assertElementsAlmostEqual(full(dij.mAlphaDose{1}), ...
        2*[1 0; 0 1],'absolute',1e-12);
    assertElementsAlmostEqual(full(dij.mAlphaDose{3}), ...
        3*[2 0; 0 2],'absolute',1e-12);
    assertElementsAlmostEqual(full(dij.mSqrtBetaDose{1}), ...
        4*[1 0; 0 1],'absolute',1e-12);
    assertElementsAlmostEqual(full(dij.mSqrtBetaDose{3}), ...
        5*[2 0; 0 2],'absolute',1e-12);
    assertElementsAlmostEqual(full(dij.customQuantity{1}), ...
        8*[1 0; 0 1],'absolute',1e-12);
    assertElementsAlmostEqual(full(dij.customQuantity{3}), ...
        9*[2 0; 0 2],'absolute',1e-12);

function test_assembler_inserts_biological_ct_metadata
    scenarioModel = fixtureTwoCtScenarioModel();
    dij1 = scenarioDij([1 0; 0 1]);
    dij2 = scenarioDij([2 0; 0 2]);
    dij1 = addBiologicalCtMetadata(dij1,[11;12],[],[21;22],[], ...
        [31;32],[],[41;42],[],[1;2],[],[3;4],[]);
    dij2 = addBiologicalCtMetadata(dij2,[],[13;14],[],[23;24], ...
        [],[33;34],[],[43;44],[],[1;2],[],[7;8]);

    dij = matRad_assembleParallelScenarioDij({dij1;dij2},[1;2], ...
        scenarioModel);

    assertEqual(dij.ax{1},[11;12]);
    assertEqual(dij.ax{2},[13;14]);
    assertEqual(dij.bx{1},[21;22]);
    assertEqual(dij.bx{2},[23;24]);
    assertEqual(dij.abx{1},[31;32]);
    assertEqual(dij.abx{2},[33;34]);
    assertEqual(dij.gamma{1},[41;42]);
    assertEqual(dij.gamma{2},[43;44]);
    assertEqual(dij.ixDose{1},[1;2]);
    assertEqual(dij.ixDose{2},[1;2]);
    assertEqual(dij.vTissueIndex{1},[3;4]);
    assertEqual(dij.vTissueIndex{2},[7;8]);

function test_assembler_accepts_identical_biological_metadata_for_shared_ct
    scenarioModel = fixtureScenarioModel();
    dij1 = scenarioDij([1 0; 0 1]);
    dij2 = scenarioDij([2 0; 0 2]);
    dij1.ax = {[1;2]};
    dij2.ax = {[1;2]};
    dij1.vTissueIndex = {[3;4]};
    dij2.vTissueIndex = {[3;4]};

    dij = matRad_assembleParallelScenarioDij({dij1;dij2},[1;2], ...
        scenarioModel);

    assertEqual(dij.ax{1},[1;2]);
    assertEqual(dij.vTissueIndex{1},[3;4]);

function test_assembler_rejects_divergent_biological_metadata_for_shared_ct
    scenarioModel = fixtureScenarioModel();
    dij1 = scenarioDij([1 0; 0 1]);
    dij2 = scenarioDij([2 0; 0 2]);
    dij1.ax = {[1;2]};
    dij2.ax = {[2;3]};

    assertExceptionThrown(@() matRad_assembleParallelScenarioDij( ...
        {dij1;dij2},[1;2],scenarioModel),'matRad:Error');

function test_assembler_rejects_biological_metadata_with_wrong_length
    scenarioModel = fixtureScenarioModel();
    dij1 = scenarioDij([1 0; 0 1]);
    dij2 = scenarioDij([2 0; 0 2]);
    dij1.ax = {[1;2]};
    dij2.ax = {[1;2;3]};

    assertExceptionThrown(@() matRad_assembleParallelScenarioDij( ...
        {dij1;dij2},[1;2],scenarioModel),'matRad:Error');

function test_assembler_rejects_biological_metadata_with_invalid_indices
    scenarioModel = fixtureScenarioModel();
    dij1 = scenarioDij([1 0; 0 1]);
    dij2 = scenarioDij([2 0; 0 2]);
    dij1.ixDose = {[1;2]};
    dij2.ixDose = {[1;3]};

    assertExceptionThrown(@() matRad_assembleParallelScenarioDij( ...
        {dij1;dij2},[1;2],scenarioModel),'matRad:Error');

function test_assembler_rejects_incompatible_grids
    scenarioModel = fixtureScenarioModel();
    dij1 = scenarioDij([1 0; 0 1]);
    dij2 = scenarioDij([2 0; 0 2]);
    dij2.doseGrid.numOfVoxels = 3;

    assertExceptionThrown(@() matRad_assembleParallelScenarioDij( ...
        {dij1;dij2},[1;2],scenarioModel),'matRad:Error');

function test_assembler_rejects_incompatible_bixel_count
    scenarioModel = fixtureScenarioModel();
    dij1 = scenarioDij([1 0; 0 1]);
    dij2 = scenarioDij([2 0; 0 2]);
    dij2.totalNumOfBixels = 3;

    assertExceptionThrown(@() matRad_assembleParallelScenarioDij( ...
        {dij1;dij2},[1;2],scenarioModel),'matRad:Error');

function test_parallel_pencil_beam_multiscen_matches_serial_or_fallback
    [ct,cst,pln,stf] = photonTestDataFixture();
    pln.multScen = rangeRandomScenario(ct);

    pln.propDoseCalc.UseParallel = false;
    dijSerial = matRad_calcDoseInfluence(ct,cst,stf,pln);

    pln.propDoseCalc.UseParallel = true;
    dijParallel = matRad_calcDoseInfluence(ct,cst,stf,pln);

    assertDijAlmostEqual(dijSerial,dijParallel);

function test_parallel_pencil_beam_multiscen_uses_parfor_when_available
    if ~parallelComputingAvailable()
        moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is unavailable.');
    end

    [ct,cst,pln,stf] = photonTestDataFixture();
    pln.multScen = rangeRandomScenario(ct);

    pln.propDoseCalc.UseParallel = false;
    dijSerial = matRad_calcDoseInfluence(ct,cst,stf,pln);

    pln.propDoseCalc.UseParallel = true;
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
    [dijParallel,useParallel] = matRad_calcDoseInfluenceParallelScenarios( ...
        ct,cst,stf,pln,engine);

    if ~useParallel
        moxunit_throw_test_skipped_exception(['Parallel scenario dij ', ...
            'calculation could not be activated in this environment.']);
    end
    assertDijAlmostEqual(dijSerial,dijParallel);

function test_parallel_pencil_beam_multiscen_engine_handle_is_worker_safe
    if ~parallelComputingAvailable()
        moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is unavailable.');
    end

    [ct,cst,pln,stf] = photonTestDataFixture();
    pln.multScen = rangeRandomScenario(ct);

    pln.propDoseCalc.UseParallel = false;
    dijSerial = matRad_calcDoseInfluence(ct,cst,stf,pln);

    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
    engine.UseParallel = true;
    engine.randomSeed = 19;
    originalRandomSeed = engine.randomSeed;
    originalScenarioIds = engine.multScen.scenarioIds();

    pln.propDoseCalc = engine;
    [dijParallel,useParallel] = matRad_calcDoseInfluenceParallelScenarios( ...
        ct,cst,stf,pln,engine);

    assertTrue(engine.UseParallel);
    assertEqual(engine.randomSeed,originalRandomSeed);
    assertEqual(engine.multScen.scenarioIds(),originalScenarioIds);
    if ~useParallel
        moxunit_throw_test_skipped_exception(['Parallel scenario dij ', ...
            'calculation could not be activated in this environment.']);
    end
    assertDijAlmostEqual(dijSerial,dijParallel);

function test_parallel_particle_bio_multiscen_matches_serial_when_available
    if ~parallelComputingAvailable()
        moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is unavailable.');
    end

    [ct,cst,pln,stf] = particleBioTestDataFixture();
    pln.multScen = rangeRandomScenario(ct);

    pln.propDoseCalc.UseParallel = false;
    dijSerial = matRad_calcDoseInfluence(ct,cst,stf,pln);

    pln.propDoseCalc.UseParallel = true;
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
    [dijParallel,useParallel] = matRad_calcDoseInfluenceParallelScenarios( ...
        ct,cst,stf,pln,engine);

    if ~useParallel
        moxunit_throw_test_skipped_exception(['Parallel particle dij ', ...
            'calculation could not be activated in this environment.']);
    end
    assertDijAlmostEqual(dijSerial,dijParallel);

function test_parallel_particle_lem_multiscen_matches_serial_when_available
    if ~parallelComputingAvailable()
        moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is unavailable.');
    end

    [ct,cst,pln,stf] = particleLemTestDataFixture();
    pln.multScen = rangeRandomScenario(ct);

    pln.propDoseCalc.UseParallel = false;
    dijSerial = matRad_calcDoseInfluence(ct,cst,stf,pln);

    pln.propDoseCalc.UseParallel = true;
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);
    [dijParallel,useParallel] = matRad_calcDoseInfluenceParallelScenarios( ...
        ct,cst,stf,pln,engine);

    if ~useParallel
        moxunit_throw_test_skipped_exception(['Parallel particle LEM dij ', ...
            'calculation could not be activated in this environment.']);
    end
    assertTrue(isfield(dijSerial,'vTissueIndex'));
    assertDijAlmostEqual(dijSerial,dijParallel);

function test_parallel_helper_falls_back_for_single_scenario
    [ct,cst,pln,stf] = photonTestDataFixture();
    pln.multScen = matRad_NominalScenario(ct);
    pln.propDoseCalc.UseParallel = true;
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

    [dij,useParallel] = matRad_calcDoseInfluenceParallelScenarios( ...
        ct,cst,stf,pln,engine);

    assertFalse(useParallel);
    assertTrue(isempty(dij));

function test_parallel_helper_falls_back_for_non_pencil_beam_engine
    pln = struct('radiationMode','brachy','machine','Generic');
    pln.propDoseCalc.UseParallel = true;
    engine = DoseEngines.matRad_TG43BrachyEngine(pln);

    [dij,useParallel] = matRad_calcDoseInfluenceParallelScenarios( ...
        [],[],[],pln,engine);

    assertFalse(useParallel);
    assertTrue(isempty(dij));

function test_parallel_helper_falls_back_without_parallel_toolbox
    if parallelComputingAvailable()
        moxunit_throw_test_skipped_exception('Parallel Computing Toolbox is available.');
    end

    [ct,cst,pln,stf] = photonTestDataFixture();
    pln.multScen = rangeRandomScenario(ct);
    pln.propDoseCalc.UseParallel = true;
    engine = DoseEngines.matRad_DoseEngineBase.getEngineFromPln(pln);

    [dij,useParallel] = matRad_calcDoseInfluenceParallelScenarios( ...
        ct,cst,stf,pln,engine);

    assertFalse(useParallel);
    assertTrue(isempty(dij));

function scenarioModel = fixtureScenarioModel()
    ct.numOfCtScen = 1;
    ct.cubeDim = [2 1 1];
    ct.resolution = struct('x',1,'y',1,'z',1);
    scenarioModel = matRad_NominalScenario(ct);
    dimensions = matRad_createScenarioComponents([1 1 1],1,1);
    scenarioValues = [0 0 0 0 0; 1 0 0 0 0];
    ctScenIds = [1;1];
    linearMask = [1 1 1; 1 3 1];
    scenMask = false(1,3,1);
    scenMask(1,1,1) = true;
    scenMask(1,3,1) = true;
    scenarioWeights = [0.5;0.5];
    scenForProb = [ctScenIds scenarioValues];
    scenarioModel.setScenarioRealizations(dimensions,scenarioValues,ctScenIds, ...
        scenarioWeights,scenarioWeights,scenForProb,linearMask,scenMask);

function scenarioModel = fixtureTwoCtScenarioModel()
    ct.numOfCtScen = 2;
    ct.cubeDim = [2 1 1];
    ct.resolution = struct('x',1,'y',1,'z',1);
    scenarioModel = matRad_NominalScenario(ct);
    dimensions = matRad_createScenarioComponents([1 1 1],1,1);
    scenarioValues = [0 0 0 0 0; 0 0 0 0 0];
    ctScenIds = [1;2];
    linearMask = [1 1 1; 2 1 1];
    scenMask = false(2,1,1);
    scenMask(1,1,1) = true;
    scenMask(2,1,1) = true;
    scenarioWeights = [0.5;0.5];
    scenForProb = [ctScenIds scenarioValues];
    scenarioModel.setScenarioRealizations(dimensions,scenarioValues,ctScenIds, ...
        scenarioWeights,scenarioWeights,scenForProb,linearMask,scenMask);

function dij = scenarioDij(matrix)
    matrix = sparse(matrix);
    dij.ctGrid = doseGrid(size(matrix,1));
    dij.doseGrid = dij.ctGrid;
    dij.numOfBeams = 1;
    dij.numOfScenarios = 1;
    dij.scenarioModel = matRad_NominalScenario();
    dij.scenarioIds = 1;
    dij.numOfRaysPerBeam = size(matrix,2);
    dij.totalNumOfBixels = size(matrix,2);
    dij.totalNumOfRays = size(matrix,2);
    dij.beamNum = ones(size(matrix,2),1);
    dij.rayNum = (1:size(matrix,2))';
    dij.bixelNum = ones(size(matrix,2),1);
    dij.minMU = zeros(size(matrix,2),1);
    dij.maxMU = inf(size(matrix,2),1);
    dij.numOfParticlesPerMU = 1e6*ones(size(matrix,2),1);
    dij.physicalDose = {matrix};

function grid = doseGrid(numVoxels)
    grid.dimensions = [numVoxels 1 1];
    grid.numOfVoxels = numVoxels;
    grid.resolution = struct('x',1,'y',1,'z',1);
    grid.x = 1;
    grid.y = 1:numVoxels;
    grid.z = 1;

function [ct,cst,pln,stf] = photonTestDataFixture()
    testDataPath = fullfile(fileparts(mfilename('fullpath')), ...
        '..','testData','photons_testData.mat');
    data = load(testDataPath,'ct','cst','pln','stf');
    ct = data.ct;
    cst = data.cst;
    pln = data.pln;
    stf = data.stf;
    pln.propDoseCalc.engine = 'SVDPB';

function [ct,cst,pln,stf] = particleBioTestDataFixture()
    testDataPath = fullfile(fileparts(mfilename('fullpath')), ...
        '..','testData','protons_testData.mat');
    data = load(testDataPath,'ct','cst','pln','stf');
    ct = data.ct;
    cst = data.cst;
    pln = data.pln;
    stf = data.stf;
    pln.propDoseCalc.engine = 'HongPB';
    pln.bioParam = matRad_bioModel(pln.radiationMode,'RBExD','MCN');

function [ct,cst,pln,stf] = particleLemTestDataFixture()
    testDataPath = fullfile(fileparts(mfilename('fullpath')), ...
        '..','testData','carbon_testData.mat');
    data = load(testDataPath,'ct','cst','pln','stf');
    ct = data.ct;
    cst = data.cst;
    pln = data.pln;
    stf = data.stf;
    pln.propDoseCalc.engine = 'HongPB';
    pln.bioParam = matRad_bioModel(pln.radiationMode,'RBExD','LEM');

function dij = addBiologicalCtMetadata(dij,ax1,ax2,bx1,bx2,abx1,abx2, ...
        gamma1,gamma2,ixDose1,ixDose2,vTissueIndex1,vTissueIndex2)
    dij.ax = {ax1;ax2};
    dij.bx = {bx1;bx2};
    dij.abx = {abx1;abx2};
    dij.gamma = {gamma1;gamma2};
    dij.ixDose = {ixDose1;ixDose2};
    dij.vTissueIndex = {vTissueIndex1;vTissueIndex2};

function scenario = rangeRandomScenario(ct)
    scenario = matRad_RandomScenarios(ct);
    scenario.nSamples = 2;
    scenario.randomSeed = 31;
    scenario.scenarioDimensionActive = {'ct','range'};

function assertDijAlmostEqual(expectedDij,actualDij)
    assertEqual(actualDij.scenarioIds,expectedDij.scenarioIds);
    assertEqual(actualDij.numOfScenarios,expectedDij.numOfScenarios);
    assertEqual(actualDij.totalNumOfBixels,expectedDij.totalNumOfBixels);
    assertEqual(actualDij.beamNum,expectedDij.beamNum);
    assertEqual(actualDij.rayNum,expectedDij.rayNum);
    assertEqual(actualDij.bixelNum,expectedDij.bixelNum);
    matrixFields = {'physicalDose','mLETDose','mAlphaDose','mSqrtBetaDose'};
    for fieldIx = 1:numel(matrixFields)
        assertDijCellFieldAlmostEqual(expectedDij,actualDij, ...
            matrixFields{fieldIx},1e-10);
    end

    bioFields = {'ax','bx','abx','gamma','ixDose','vTissueIndex'};
    for fieldIx = 1:numel(bioFields)
        assertDijCellFieldAlmostEqual(expectedDij,actualDij, ...
            bioFields{fieldIx},1e-12);
    end

function assertDijCellFieldAlmostEqual(expectedDij,actualDij,fieldName,tolerance)
    if ~isfield(expectedDij,fieldName)
        assertFalse(isfield(actualDij,fieldName));
        return;
    end

    assertTrue(isfield(actualDij,fieldName));
    assertEqual(size(actualDij.(fieldName)),size(expectedDij.(fieldName)));

    for i = 1:numel(expectedDij.(fieldName))
        if isempty(expectedDij.(fieldName){i})
            assertTrue(isempty(actualDij.(fieldName){i}));
        else
            assertElementsAlmostEqual(full(actualDij.(fieldName){i}), ...
                full(expectedDij.(fieldName){i}),'absolute',tolerance);
        end
    end

function available = parallelComputingAvailable()
    available = false;
    if exist('parpool','file') ~= 2 || exist('gcp','file') ~= 2
        return;
    end

    try
        [available,~] = license('checkout','Distrib_Computing_Toolbox');
    catch
        available = false;
    end

    if isempty(available)
        available = false;
    end
    available = logical(available);
