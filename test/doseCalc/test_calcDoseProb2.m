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
        full(dij.physicalDose{1}),'absolute',1e-12);
    assertTrue(any(abs(full(dijProb2Context.physicalDose{1}(:)) - ...
        expected(:)) > 1e-12));
    assertEqual(dijProb2Context.numOfScenarios,1);
    assertTrue(isa(dijProb2Context.scenarioModel,'matRad_NominalScenario'));
    assertEqual(dijProb2Context.probabilisticQuantity,'physicalDose');
    assertEqual(dijProb2Context.probabilisticQuantityField,'physicalDose');

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

function test_prob2_response_is_consistent_with_selected_vois
    [ct,cst,pln,dij,cfg] = partialSelectionFixture();

    [plnOut,~] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);
    dijProb2 = plnOut.propOpt.dij_prob2;

    expected = [2.5 0; 0 3.5; 1.75 1; 0 0];
    assertElementsAlmostEqual(full(dijProb2.expected),expected,'absolute',1e-12);
    assertEqual(dijProb2.voiSubIx{1},[1;2]);
    assertEqual(dijProb2.voiSubIx{2},3);
    assertTrue(isempty(dijProb2.voiSubIx{3}));
    assertElementsAlmostEqual(full(dijProb2.Omega{1}),diag([0.75 0.75]), ...
        'absolute',1e-12);
    assertElementsAlmostEqual(full(dijProb2.Omega{2}), ...
        [0.1875 0; 0 0],'absolute',1e-12);
    assertTrue(isempty(dijProb2.Omega{3}));

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

function test_prob2_streaming_inmemory_matches_existing
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    [plnLegacy,~] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);

    [plnStream,dijStreamContext] = matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg);

    assertProb2AlmostEqual(plnLegacy.propOpt.dij_prob2, ...
        plnStream.propOpt.dij_prob2);
    assertEqual(plnStream.propOpt.dij_prob2.precomputeMode,'streaming');
    assertEqual(plnStream.propOpt.dij_prob2.secondPassStrategy,'disk');
    assertEqual(dijStreamContext.numOfScenarios,1);
    assertElementsAlmostEqual(full(dijStreamContext.physicalDose{1}), ...
        full(dij.physicalDose{1}),'absolute',1e-12);
    assertTrue(any(abs(full(dijStreamContext.physicalDose{1}(:)) - ...
        full(plnLegacy.propOpt.dij_prob2.expected(:))) > 1e-12));

function test_prob2_streaming_partial_selection_matches_inmemory
    [ct,cst,pln,dij,cfg] = partialSelectionFixture();
    cfg.BatchSize = 1;
    [plnLegacy,~] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);

    [plnStreamDisk,~] = matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg);

    assertProb2AlmostEqual(plnLegacy.propOpt.dij_prob2, ...
        plnStreamDisk.propOpt.dij_prob2);
    assertEqual(plnStreamDisk.propOpt.dij_prob2.secondPassStrategy,'disk');

    cfg.SecondPassStrategy = 'recompute';
    [plnStreamRecompute,~] = matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg);

    assertProb2AlmostEqual(plnLegacy.propOpt.dij_prob2, ...
        plnStreamRecompute.propOpt.dij_prob2);
    assertEqual(plnStreamRecompute.propOpt.dij_prob2.secondPassStrategy, ...
        'recompute');

function test_prob2_streaming_accepts_dij_without_cfg
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    [plnLegacy,~] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);

    [plnStream,~] = matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij);

    assertProb2AlmostEqual(plnLegacy.propOpt.dij_prob2, ...
        plnStream.propOpt.dij_prob2);

function test_prob2_streaming_rejects_duplicate_precomputed_dij_inputs
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.PrecomputedDij = dij;

    assertExceptionThrown(@() matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_prob2_streaming_recomputes_scenario_dij_without_dij_argument
    [ct,cst,pln,stf] = photonTestDataFixture();
    cfg.SecondPassStrategy = 'disk';
    cfg.KeepCache = false;
    cfg.BatchSize = 10000;
    cfg.targetStructSel = 1;

    [plnStream,dijStreamContext] = matRad_calcDoseProb2Streaming(ct,cst,stf,pln,cfg);
    dijProb2 = plnStream.propOpt.dij_prob2;

    assertEqual(size(dijProb2.expected,1),dijStreamContext.doseGrid.numOfVoxels);
    assertEqual(size(dijProb2.expected,2),dijStreamContext.totalNumOfBixels);
    assertTrue(nnz(dijProb2.expected) > 0);
    assertFalse(isfield(dijProb2,'cacheDir'));
    assertStreamingSizeDisk(dijProb2);
    assertEqual(dijStreamContext.numOfScenarios,1);

function test_prob2_streaming_recompute_matches_existing
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.BatchSize = 1;
    [plnLegacy,~] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);

    cfg.SecondPassStrategy = 'recompute';
    [plnStream,~] = matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg);

    assertProb2AlmostEqual(plnLegacy.propOpt.dij_prob2, ...
        plnStream.propOpt.dij_prob2);
    assertEqual(plnStream.propOpt.dij_prob2.secondPassStrategy,'recompute');
    assertStreamingSizeRecompute(plnStream.propOpt.dij_prob2);

function test_prob2_streaming_disk_cache_keeps_distinct_hash_folders
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cacheRoot = tempname();
    cleanup = onCleanup(@() deleteDirIfExists(cacheRoot));
    cfg.CacheRoot = cacheRoot;
    cfg.KeepCache = true;

    [plnStream1,~] = matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg);
    [plnStream2,~] = matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg);

    prob2_1 = plnStream1.propOpt.dij_prob2;
    prob2_2 = plnStream2.propOpt.dij_prob2;
    assertTrue(isfield(prob2_1,'cacheDir'));
    assertTrue(isfield(prob2_2,'cacheDir'));
    assertStreamingSizeDisk(prob2_1);
    assertStreamingSizeDisk(prob2_2);
    assertFalse(strcmp(prob2_1.cacheDir,prob2_2.cacheDir));
    assertEqual(exist(prob2_1.cacheDir,'dir'),7);
    assertEqual(exist(prob2_2.cacheDir,'dir'),7);
    assertEqual(exist(fullfile(prob2_1.cacheDir,'metadata.mat'),'file'),2);
    metadataData = load(fullfile(prob2_1.cacheDir,'metadata.mat'));
    assertTrue(isfield(metadataData,'metadata'));
    metadata = metadataData.metadata;
    assertEqual(metadata.calculationName,'PROB2');
    assertEqual(metadata.quantity,'physicalDose');
    assertEqual(metadata.refScen,1);
    assertEqual(metadata.scenarioDijIx,[1;2]);
    assertEqual(metadata.scenarioCtScenIds,[1;1]);
    assertElementsAlmostEqual(metadata.scenarioWeights,[0.25;0.75], ...
        'absolute',1e-12);
    assertTrue(~isempty(dir(fullfile(prob2_1.cacheDir,'scenario_*_voi_*_block_*.mat'))));

function test_prob2_streaming_disk_cache_cleans_hash_folder_by_default
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cacheRoot = tempname();
    cleanup = onCleanup(@() deleteDirIfExists(cacheRoot));
    cfg.CacheRoot = cacheRoot;
    cfg.KeepCache = false;

    [plnStream,~] = matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg);

    dijProb2 = plnStream.propOpt.dij_prob2;
    assertFalse(isfield(dijProb2,'cacheDir'));
    assertStreamingSizeDisk(dijProb2);
    assertEqual(numel(listCacheRunDirs(cacheRoot)),0);

function test_prob2_streaming_multict_mapping_matches_existing
    [ct,cst,pln,dij,cfg] = multiCtFixture();
    pln.propOpt.scen4D = 'all';
    [plnLegacy,~] = matRad_calcDoseProb2(ct,cst,[],pln,dij,cfg);

    [plnStream,~] = matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg);

    assertProb2AlmostEqual(plnLegacy.propOpt.dij_prob2, ...
        plnStream.propOpt.dij_prob2);
    assertEqual(plnStream.propOpt.dij_prob2.scenarioCtScenIds,[1;2]);

function test_prob2_streaming_rejects_invalid_second_pass_strategy
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.SecondPassStrategy = 'memory';

    assertExceptionThrown(@() matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_prob2_streaming_rejects_invalid_cache_root
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.CacheRoot = 1;

    assertExceptionThrown(@() matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_prob2_streaming_rejects_cache_root_file
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.CacheRoot = tempname();
    fid = fopen(cfg.CacheRoot,'w');
    fwrite(fid,'not a folder');
    fclose(fid);
    cleanup = onCleanup(@() deleteFileIfExists(cfg.CacheRoot));

    assertExceptionThrown(@() matRad_calcDoseProb2Streaming(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

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

function [ct,cst,pln,dij,cfg] = partialSelectionFixture()
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cst{2,4}{1} = 3;
    cst(3,1:6) = cell(1,6);
    cst = addStructure(cst,3,'Spare OAR','OAR',4);
    cfg.OARStructSel = 2;

function [ct,cst,pln,dij,cfg] = multiCtFixture()
    dim = [2 3 1];
    [xGrid,~,~] = meshgrid(1:dim(2),1:dim(1),1:dim(3));
    sourceRows = xGrid(:);

    ct.numOfCtScen = 2;
    ct.cubeDim = dim;
    ct.resolution = struct('x',1,'y',1,'z',1);
    ct.x = 1:dim(2);
    ct.y = 1:dim(1);
    ct.z = 1:dim(3);
    ct.refScen = 1;
    ct.dvfMetadata.dvfType = 'pull';
    ct.dvfMetadata.dvfUnits = 'voxel';
    ct.dvfMetadata.refScen = 1;
    ct.dvfMetadata.referenceCtScen = 1;
    ct.dvf = cell(2,1);
    ct.dvf{1} = zeros([3 dim]);
    ct.dvf{2} = zeros([3 dim]);
    ct.dvf{2}(1,:,:,:) = 1;

    cst = cell(2,6);
    cst = addStructure(cst,1,'PTV','TARGET',[1;2;3]);
    cst = addStructure(cst,2,'OAR','OAR',[4;5;6]);
    cst{1,4} = cell(1,2);
    cst{2,4} = cell(1,2);
    cst{1,4}{1} = [1;2;3];
    cst{2,4}{1} = [4;5;6];

    pln.bioParam.quantityOpt = 'physicalDose';
    pln.multScen = fixtureScenarioModel(ct,[1;2],zeros(2,5), ...
        [1 1 1; 2 1 1],true(2,1,1),[0.5;0.5]);

    dij = baseDij(dim,1);
    dij.physicalDose = cell(2,1,1);
    dij.physicalDose{1} = sparse(sourceRows);
    dij.physicalDose{2} = sparse(sourceRows);
    cfg.refScen = 1;
    cfg.BatchSize = 2;

function [ct,cst,pln,stf] = photonTestDataFixture()
    testDataPath = fullfile(fileparts(mfilename('fullpath')), ...
        '..','testData','photons_testData.mat');
    data = load(testDataPath,'ct','cst','pln','stf');
    ct = data.ct;
    cst = data.cst;
    pln = data.pln;
    stf = data.stf;
    pln.propDoseCalc.engine = 'SVDPB';
    if ~isfield(pln,'multScen') || isempty(pln.multScen)
        pln.multScen = matRad_NominalScenario();
    end

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

function assertProb2AlmostEqual(expectedProb2,actualProb2)
    assertElementsAlmostEqual(full(actualProb2.expected), ...
        full(expectedProb2.expected),'absolute',1e-12);
    assertEqual(size(actualProb2.Omega),size(expectedProb2.Omega));
    for i = 1:numel(expectedProb2.Omega)
        if isempty(expectedProb2.Omega{i})
            assertTrue(isempty(actualProb2.Omega{i}));
        else
            assertElementsAlmostEqual(full(actualProb2.Omega{i}), ...
                full(expectedProb2.Omega{i}),'absolute',1e-12);
        end
    end
    assertEqual(actualProb2.voiSubIx,expectedProb2.voiSubIx);
    assertEqual(actualProb2.quantity,expectedProb2.quantity);
    assertEqual(actualProb2.quantityField,expectedProb2.quantityField);
    assertElementsAlmostEqual(actualProb2.scenarioWeights, ...
        expectedProb2.scenarioWeights,'absolute',1e-12);

function assertStreamingSizeTotal(dij)
    assertTrue(isfield(dij,'streamingSize'));
    sizeData = dij.streamingSize;
    assertTrue(sizeData.compactBytes > 0);
    assertTrue(sizeData.auxiliaryPeakBytes >= 0);
    assertElementsAlmostEqual(sizeData.totalPrecomputingBytes, ...
        sizeData.compactBytes + sizeData.auxiliaryPeakBytes, ...
        'absolute',1e-12);

function assertStreamingSizeDisk(dij)
    assertStreamingSizeTotal(dij);
    sizeData = dij.streamingSize;
    assertTrue(sizeData.diskCachePeakBytes > 0);
    assertEqual(sizeData.auxiliaryPeakBytes, ...
        sizeData.diskCachePeakBytes);
    assertEqual(sizeData.auxiliaryPeakKind,'diskCache');
    assertEqual(sizeData.secondPassStrategy,'disk');

function assertStreamingSizeRecompute(dij)
    assertStreamingSizeTotal(dij);
    sizeData = dij.streamingSize;
    assertTrue(sizeData.memoryTemporaryPeakBytes > 0);
    assertEqual(sizeData.auxiliaryPeakBytes, ...
        sizeData.memoryTemporaryPeakBytes);
    assertEqual(sizeData.auxiliaryPeakKind,'memoryTemporary');
    assertEqual(sizeData.secondPassStrategy,'recompute');

function dirs = listCacheRunDirs(cacheRoot)
    dirs = {};
    if exist(cacheRoot,'dir') ~= 7
        return;
    end
    listing = dir(cacheRoot);
    isRunDir = [listing.isdir] & ~ismember({listing.name},{'.','..'});
    dirs = {listing(isRunDir).name};

function deleteDirIfExists(path)
    if exist(path,'dir') == 7
        rmdir(path,'s');
    end

function deleteFileIfExists(path)
    if exist(path,'file') == 2
        delete(path);
    end
