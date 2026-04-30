function test_suite = test_calcDoseInterval

test_functions = localfunctions();

initTestSuite;

function test_interval2_physical_dose_center_and_radius
    [ct,cst,pln,dij,cfg] = singleCtFixture();

    [~,~,~,plnOut,dijInterval] = matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(dijInterval.center(1:4,:)), ...
        [2.5 0; 0 3.5; 1.75 1; 0 4],'absolute',1e-12);
    assertElementsAlmostEqual(full(dijInterval.radius),diag([7 13]),'absolute',1e-12);
    assertEqual(dijInterval.quantity,'physicalDose');
    assertEqual(dijInterval.quantityField,'physicalDose');
    assertEqual(dijInterval.refScen,1);
    assertElementsAlmostEqual(dijInterval.scenarioWeights,[0.25;0.75],'absolute',1e-12);
    assertTrue(isfield(plnOut.propOpt,'dij_interval'));

function test_interval2_explicit_linear_quantity
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    dij.mAlphaDose = cell(size(dij.physicalDose));
    dij.mAlphaDose{1} = 2*dij.physicalDose{1};
    dij.mAlphaDose{2} = 2*dij.physicalDose{2};
    cfg.Quantity = 'mAlphaDose';

    [~,~,~,~,dijInterval] = matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(dijInterval.center(1:2,:)), ...
        2*[2.5 0; 0 3.5],'absolute',1e-12);
    assertElementsAlmostEqual(full(dijInterval.radius),4*diag([7 13]),'absolute',1e-12);
    assertEqual(dijInterval.quantityField,'mAlphaDose');

function test_interval2_const_rbe_scales_quantity
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    dij.RBE = 1.1;
    cfg.Quantity = 'RBExD';

    [~,~,~,~,dijInterval] = matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(dijInterval.center(1:2,:)), ...
        1.1*[2.5 0; 0 3.5],'absolute',1e-12);
    assertElementsAlmostEqual(full(dijInterval.radius),1.21*diag([7 13]),'absolute',1e-12);
    assertEqual(dijInterval.quantity,'RBExD');
    assertEqual(dijInterval.quantityField,'physicalDose');

function test_interval_rejects_nonlinear_quantity_without_linear_field
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.Quantity = 'effect';

    assertExceptionThrown(@() matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_interval_rejects_invalid_progress_level
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.ProgressLevel = 'verbose';

    assertExceptionThrown(@() matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_interval2_batch_size_does_not_change_result
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.BatchSize = 1;
    [~,~,~,~,dijIntervalBatch] = matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg);

    cfg.BatchSize = 99;
    [~,~,~,~,dijIntervalFull] = matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(dijIntervalBatch.center), ...
        full(dijIntervalFull.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(dijIntervalBatch.radius), ...
        full(dijIntervalFull.radius),'absolute',1e-12);

function test_interval3_oar_covariance_svd
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.KMode = 'static';
    cfg.KMax = 2;

    [~,~,~,~,dijInterval] = matRad_calcDoseInterval3(ct,cst,[],pln,dij,cfg);

    covRow3 = full(dijInterval.U{1}*dijInterval.S{1}*dijInterval.V{1}');
    covRow4 = full(dijInterval.U{2}*dijInterval.S{2}*dijInterval.V{2}');
    assertElementsAlmostEqual(covRow3,[0.1875 0; 0 0],'absolute',1e-10);
    assertElementsAlmostEqual(covRow4,[0 0; 0 3],'absolute',1e-10);
    assertEqual(dijInterval.OARSubIx,[3;4]);

function test_interval3_oar_svd_accepts_sufficient_memory_limit
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.KMode = 'static';
    cfg.KMax = 2;
    [~,~,~,~,dijIntervalDefault] = matRad_calcDoseInterval3(ct,cst,[],pln,dij,cfg);

    cfg.MemoryLimitMB = 1;
    [~,~,~,~,dijIntervalLimited] = matRad_calcDoseInterval3(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(dijIntervalLimited.center), ...
        full(dijIntervalDefault.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(dijIntervalLimited.U{1}*dijIntervalLimited.S{1}*dijIntervalLimited.V{1}'), ...
        full(dijIntervalDefault.U{1}*dijIntervalDefault.S{1}*dijIntervalDefault.V{1}'),'absolute',1e-12);

function test_interval3_oar_svd_rejects_low_memory_limit
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.MemoryLimitMB = 1e-6;

    assertExceptionThrown(@() matRad_calcDoseInterval3(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_interval2_low_memory_limit_does_not_apply_oar_svd_guard
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.MemoryLimitMB = 1e-6;

    [~,~,~,~,dijInterval] = matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(dijInterval.center(1:4,:)), ...
        [2.5 0; 0 3.5; 1.75 1; 0 4],'absolute',1e-12);

function test_interval3_batch_size_does_not_change_result_when_memory_allows
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.KMode = 'static';
    cfg.KMax = 2;
    cfg.MemoryLimitMB = 1;
    cfg.BatchSize = 1;
    [~,~,~,~,dijIntervalBatch] = matRad_calcDoseInterval3(ct,cst,[],pln,dij,cfg);

    cfg.BatchSize = 99;
    [~,~,~,~,dijIntervalFull] = matRad_calcDoseInterval3(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(dijIntervalBatch.center), ...
        full(dijIntervalFull.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(dijIntervalBatch.U{2}*dijIntervalBatch.S{2}*dijIntervalBatch.V{2}'), ...
        full(dijIntervalFull.U{2}*dijIntervalFull.S{2}*dijIntervalFull.V{2}'),'absolute',1e-12);

function test_interval3_zero_oar_covariance_is_valid_zero_rank
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    dij.physicalDose{2}(3,:) = dij.physicalDose{1}(3,:);
    cfg.OARStructSel = 'OAR';
    cfg.targetStructSel = 'PTV';

    [~,~,~,~,dijInterval] = matRad_calcDoseInterval3(ct,cst,[],pln,dij,cfg);

    assertEqual(dijInterval.k(1),0);
    assertEqual(size(dijInterval.U{1}),[2 0]);
    assertEqual(size(dijInterval.S{1}),[0 0]);
    assertEqual(size(dijInterval.V{1}),[2 0]);

function test_interval2_multict_maps_to_reference_scenario
    [ct,cst,pln,dij,cfg,expectedCenter] = multiCtFixture(1);

    [~,~,~,~,dijInterval] = matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(dijInterval.center(:,1)),expectedCenter,'absolute',1e-12);
    assertEqual(dijInterval.refScen,1);
    assertEqual(dijInterval.scenarioCtScen,[1;2]);

function test_interval2_multict_supports_nonfirst_reference_scenario
    [ct,cst,pln,dij,cfg,expectedCenter] = multiCtFixture(2);

    [~,~,~,~,dijInterval] = matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(dijInterval.center(:,1)),expectedCenter,'absolute',1e-12);
    assertEqual(dijInterval.refScen,2);
    assertEqual(dijInterval.scenarioCtScen,[2;3]);

function test_interval2_multict_requires_pull_dvf
    [ct,cst,pln,dij,cfg] = multiCtFixture(1);
    ct.dvf = {};

    assertExceptionThrown(@() matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_interval2_multict_uses_dij_ct_grid_when_ct_axes_are_missing
    [ct,cst,pln,dij,cfg,expectedCenter] = multiCtResamplingFixture();

    [~,~,~,~,dijInterval] = matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(dijInterval.center(:,1)),expectedCenter,'absolute',1e-12);

function [ct,cst,pln,dij,cfg] = singleCtFixture()
    ct.numOfCtScen = 1;
    ct.cubeDim = [2 2 1];
    ct.resolution = struct('x',1,'y',1,'z',1);
    ct.refScen = 1;

    cst = cell(2,6);
    cst = addStructure(cst,1,'PTV','TARGET',[1;2]);
    cst = addStructure(cst,2,'OAR','OAR',[3;4]);

    pln.bioParam.quantityOpt = 'physicalDose';
    pln.multScen.scenMask = true(1,2,1);
    pln.multScen.linearMask = [1 1 1; 1 2 1];
    pln.multScen.scenWeight = [0.25;0.75];
    pln.multScen.scenProb = [0.25;0.75];

    dij = baseDij(ct.cubeDim,2);
    dij.physicalDose = cell(1,2,1);
    dij.physicalDose{1} = sparse([1 0; 0 2; 1 1; 0 1]);
    dij.physicalDose{2} = sparse([3 0; 0 4; 2 1; 0 5]);

    cfg.CalculateReferenceDij = false;
    cfg.BatchSize = 2;

function [ct,cst,pln,dij,cfg,expectedCenter] = multiCtFixture(refScen)
    dim = [2 3 1];
    [xGrid,~,~] = meshgrid(1:dim(2),1:dim(1),1:dim(3));
    sourceRows = xGrid(:);
    mappedRows = zeros(dim);
    mappedRows(:,2:end,:) = xGrid(:,1:end-1,:);
    mappedRows = mappedRows(:);

    ct.numOfCtScen = refScen + 1;
    ct.cubeDim = dim;
    ct.resolution = struct('x',1,'y',1,'z',1);
    ct.x = 1:dim(2);
    ct.y = 1:dim(1);
    ct.z = 1:dim(3);
    ct.refScen = refScen;
    ct.dvfMetadata.dvfType = 'pull';
    ct.dvfMetadata.dvfUnits = 'voxel';
    ct.dvfMetadata.refScen = refScen;
    ct.dvfMetadata.referenceCtScen = refScen;
    ct.dvf = cell(ct.numOfCtScen,1);
    ct.dvf{refScen} = zeros([3 dim]);
    ct.dvf{refScen + 1} = zeros([3 dim]);
    ct.dvf{refScen + 1}(1,:,:,:) = 1;

    cst = cell(1,6);
    cst = addStructure(cst,1,'PTV','TARGET',(1:prod(dim))');
    cst{1,4} = cell(1,ct.numOfCtScen);
    cst{1,4}{refScen} = (1:prod(dim))';

    pln.bioParam.quantityOpt = 'physicalDose';
    pln.multScen.scenMask = false(ct.numOfCtScen,1,1);
    pln.multScen.scenMask(refScen) = true;
    pln.multScen.scenMask(refScen + 1) = true;
    pln.multScen.linearMask = [refScen 1 1; refScen + 1 1 1];
    pln.multScen.scenWeight = [0.5;0.5];
    pln.multScen.scenProb = [0.5;0.5];

    dij = baseDij(dim,1);
    dij.physicalDose = cell(ct.numOfCtScen,1,1);
    dij.physicalDose{refScen} = sparse(sourceRows);
    dij.physicalDose{refScen + 1} = sparse(sourceRows);

    cfg.CalculateReferenceDij = false;
    cfg.refScen = refScen;
    cfg.BatchSize = 2;

    expectedCenter = 0.5*sourceRows + 0.5*mappedRows;

function [ct,cst,pln,dij,cfg,expectedCenter] = multiCtResamplingFixture()
    ctDim = [2 4 2];
    doseDim = [2 2 2];

    ct.numOfCtScen = 2;
    ct.cubeDim = ctDim;
    ct.resolution = struct('x',1,'y',1,'z',1);
    ct.refScen = 1;
    ct.dvfMetadata.dvfType = 'pull';
    ct.dvfMetadata.dvfUnits = 'mm';
    ct.dvfMetadata.refScen = 1;
    ct.dvfMetadata.referenceCtScen = 1;
    ct.dvf = cell(2,1);
    ct.dvf{1} = zeros([3 ctDim]);
    ct.dvf{2} = zeros([3 ctDim]);

    cst = cell(1,6);
    cst = addStructure(cst,1,'PTV','TARGET',(1:prod(ctDim))');

    pln.bioParam.quantityOpt = 'physicalDose';
    pln.multScen.scenMask = true(2,1,1);
    pln.multScen.linearMask = [1 1 1; 2 1 1];
    pln.multScen.scenWeight = [0.5;0.5];
    pln.multScen.scenProb = [0.5;0.5];

    dij = baseDij(doseDim,1);
    dij.ctGrid.dimensions = ctDim;
    dij.ctGrid.numOfVoxels = prod(ctDim);
    dij.ctGrid.resolution = ct.resolution;
    dij.ctGrid.x = 1:ctDim(2);
    dij.ctGrid.y = 1:ctDim(1);
    dij.ctGrid.z = 1:ctDim(3);
    dij.doseGrid.resolution = struct('x',2,'y',1,'z',1);
    dij.doseGrid.x = 1:2:ctDim(2);
    dij.doseGrid.y = 1:ctDim(1);
    dij.doseGrid.z = 1:ctDim(3);
    dij.physicalDose = cell(2,1,1);
    dij.physicalDose{1} = sparse((1:prod(doseDim))');
    dij.physicalDose{2} = sparse((prod(doseDim) + 1:2*prod(doseDim))');

    cfg.CalculateReferenceDij = false;
    cfg.refScen = 1;
    cfg.BatchSize = 2;
    expectedCenter = 0.5*((1:prod(doseDim))' + (prod(doseDim) + 1:2*prod(doseDim))');

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
