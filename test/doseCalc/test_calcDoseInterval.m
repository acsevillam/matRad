function test_suite = test_calcDoseInterval

test_functions = localfunctions();

initTestSuite;

function test_interval2_physical_dose_center_and_radius
    [ct,cst,pln,dij,cfg] = singleCtFixture();

    [plnOut,dijIntervalContext] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijInterval.center(1:4,:)), ...
        [2.5 0; 0 3.5; 1.75 1; 0 4],'absolute',1e-12);
    assertElementsAlmostEqual(full(dijInterval.radius),diag([7 13]),'absolute',1e-12);
    assertEqual(dijInterval.radiusMode,'std');
    assertEqual(dijInterval.intervalMode,'INTERVAL2');
    assertEqual(dijInterval.targetSubIx,[1;2]);
    assertEqual(dijInterval.OARSubIx,[3;4]);
    assertEqual(dijInterval.quantity,'physicalDose');
    assertEqual(dijInterval.quantityField,'physicalDose');
    assertEqual(dijInterval.optimizationQuantity,'physicalDose');
    assertEqual(dijInterval.refScen,1);
    assertElementsAlmostEqual(dijInterval.scenarioWeights,[0.25;0.75],'absolute',1e-12);
    assertFalse(isfield(dijInterval,'OARRadiusFactor'));
    assertTrue(isfield(plnOut.propOpt,'dij_interval'));
    assertEqual(dijIntervalContext.totalNumOfBixels,2);
    assertEqual(dijIntervalContext.numOfScenarios,1);
    assertEqual(dijIntervalContext.beamNum,ones(2,1));
    assertTrue(isa(dijIntervalContext.scenarioModel,'matRad_NominalScenario'));
    assertEqual(dijIntervalContext.scenarioModel.numScenarios(),1);
    assertEqual(dijIntervalContext.scenarioModel.getDijScenarioIndex(1),1);
    assertEqual(dijIntervalContext.scenarioModel.getCtScenario(1),1);
    assertEqual(plnOut.multScen.numScenarios(),1);
    assertEqual(plnOut.multScen.getDijScenarioIndex(1),1);
    assertElementsAlmostEqual(full(dijIntervalContext.physicalDose{1}), ...
        full(dij.physicalDose{1}),'absolute',1e-12);
    assertTrue(any(abs(full(dijIntervalContext.physicalDose{1}(:)) - ...
        full(dijInterval.center(:))) > 1e-12));

function test_interval_generic_entrypoint_uses_cfg_interval_mode
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    [plnLegacy,~] = calcInterval2(ct,cst,[],pln,dij,cfg);

    cfg.IntervalMode = 'INTERVAL2';
    [plnGeneric,dijGenericContext] = matRad_calcDoseInterval(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(plnGeneric.propOpt.dij_interval.center), ...
        full(plnLegacy.propOpt.dij_interval.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(plnGeneric.propOpt.dij_interval.radius), ...
        full(plnLegacy.propOpt.dij_interval.radius),'absolute',1e-12);
    assertEqual(plnGeneric.propOpt.dij_interval.intervalMode,'INTERVAL2');
    assertEqual(dijGenericContext.numOfScenarios,1);

function test_interval2_explicit_linear_quantity
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    dij.mAlphaDose = cell(size(dij.physicalDose));
    dij.mAlphaDose{1} = 2*dij.physicalDose{1};
    dij.mAlphaDose{2} = 2*dij.physicalDose{2};
    cfg.Quantity = 'mAlphaDose';

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijInterval.center(1:2,:)), ...
        2*[2.5 0; 0 3.5],'absolute',1e-12);
    assertElementsAlmostEqual(full(dijInterval.radius),4*diag([7 13]),'absolute',1e-12);
    assertEqual(dijInterval.quantityField,'mAlphaDose');

function test_interval2_const_rbe_scales_quantity
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    dij.RBE = 1.1;
    cfg.Quantity = 'RBExD';

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

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

function test_interval_rejects_invalid_radius_mode
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.RadiusMode = 'range';

    assertExceptionThrown(@() matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_interval2_batch_size_does_not_change_result
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.BatchSize = 1;
    [plnBatch,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijIntervalBatch = plnBatch.propOpt.dij_interval;

    cfg.BatchSize = 99;
    [plnFull,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijIntervalFull = plnFull.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijIntervalBatch.center), ...
        full(dijIntervalFull.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(dijIntervalBatch.radius), ...
        full(dijIntervalFull.radius),'absolute',1e-12);

function test_interval2_response_is_consistent_with_selected_vois
    [ct,cst,pln,dij,cfg] = partialSelectionFixture();

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    expectedCenter = [2.5 0; 0 3.5; 1.75 1; 0 0];
    assertElementsAlmostEqual(full(dijInterval.center),expectedCenter, ...
        'absolute',1e-12);
    assertElementsAlmostEqual(full(dijInterval.radius),diag([7 13]), ...
        'absolute',1e-12);
    assertEqual(dijInterval.targetSubIx,[1;2]);
    assertEqual(dijInterval.OARSubIx,3);

function test_interval_rows_follow_overlap_priorities
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cst{1,5}.Priority = 1;
    cst{2,5}.Priority = 2;
    cst{1,6} = {struct()};
    cst{2,4}{1} = [2;3;4];
    cst{2,6} = {struct()};

    [plnOut,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertEqual(dijInterval.targetSubIx,[1;2]);
    assertEqual(dijInterval.OARSubIx,[3;4]);

function test_interval2_explicit_std_matches_default_radius_mode
    [ct,cst,pln,dij,cfg] = singleCtFixture();

    [plnDefault,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    cfg.RadiusMode = 'std';
    [plnStd,~] = calcInterval2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(plnStd.propOpt.dij_interval.center), ...
        full(plnDefault.propOpt.dij_interval.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(plnStd.propOpt.dij_interval.radius), ...
        full(plnDefault.propOpt.dij_interval.radius),'absolute',1e-12);
    assertEqual(plnStd.propOpt.dij_interval.radiusMode,'std');

function test_interval2_extreme_target_radius
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.RadiusMode = 'extreme';

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertEqual(dijInterval.radiusMode,'extreme');
    assertElementsAlmostEqual(full(dijInterval.center(1:4,:)), ...
        [2.5 0; 0 3.5; 1.75 1; 0 4],'absolute',1e-12);
    assertElementsAlmostEqual(full(dijInterval.radius), ...
        diag([8.5 14.5]),'absolute',1e-12);

function test_interval2_collect_timing_reports_target_and_oar_components
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.CollectTiming = true;

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;
    timing = dijInterval.timing;

    assertEqual(timing.intervalMode,'INTERVAL2');
    assertEqual(timing.numTargetVoxels,2);
    assertEqual(timing.numOarVoxels,2);
    assertEqual(timing.numScenarios,2);
    assertEqual(timing.numBixels,2);
    assertTimingStageIsValid(timing.target);
    assertTimingStageIsValid(timing.oar);
    assertTrue(timing.target.radiusMultiplySeconds >= 0);
    assertEqual(timing.target.factorSeconds,0);
    assertEqual(timing.oar.radiusMultiplySeconds,0);
    assertEqual(timing.oar.factorSeconds,0);

function test_interval3_response_is_consistent_with_selected_vois
    [ct,cst,pln,dij,cfg] = partialSelectionFixture();
    cfg.KMode = 'static';
    cfg.KMax = 2;

    [plnOut,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    expectedCenter = [2.5 0; 0 3.5; 1.75 1; 0 0];
    assertElementsAlmostEqual(full(dijInterval.center),expectedCenter, ...
        'absolute',1e-12);
    assertElementsAlmostEqual(full(dijInterval.radius),diag([7 13]), ...
        'absolute',1e-12);
    assertEqual(dijInterval.targetSubIx,[1;2]);
    assertEqual(dijInterval.OARSubIx,3);
    assertEqual(numel(dijInterval.OARRadiusFactor),1);
    assertElementsAlmostEqual(full(reconstructOARRadiusCovariance( ...
        dijInterval,1)),[0.1875 0; 0 0],'absolute',1e-10);
    assertEqual(dijInterval.OARRadiusRank,1);

function test_interval3_oar_covariance_factor
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.KMode = 'static';
    cfg.KMax = 2;

    [plnOut,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    covRow3 = full(reconstructOARRadiusCovariance(dijInterval,1));
    covRow4 = full(reconstructOARRadiusCovariance(dijInterval,2));
    assertElementsAlmostEqual(covRow3,[0.1875 0; 0 0],'absolute',1e-10);
    assertElementsAlmostEqual(covRow4,[0 0; 0 3],'absolute',1e-10);
    assertEqual(dijInterval.radiusMode,'std');
    assertEqual(dijInterval.OARSubIx,[3;4]);
    assertEqual(dijInterval.OARRadiusRank,[1;1]);
    assertFalse(isfield(dijInterval,'U'));
    assertFalse(isfield(dijInterval,'S'));
    assertFalse(isfield(dijInterval,'V'));
    assertFalse(isfield(dijInterval,'k'));

function test_interval3_extreme_oar_radius_factor
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.RadiusMode = 'extreme';
    cfg.KMode = 'static';
    cfg.KMax = 2;

    [plnOut,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    covRow3 = full(reconstructOARRadiusCovariance(dijInterval,1));
    covRow4 = full(reconstructOARRadiusCovariance(dijInterval,2));
    assertEqual(dijInterval.radiusMode,'extreme');
    assertElementsAlmostEqual(full(dijInterval.radius), ...
        diag([8.5 14.5]),'absolute',1e-12);
    assertElementsAlmostEqual(covRow3,[0.5625 0; 0 0],'absolute',1e-10);
    assertElementsAlmostEqual(covRow4,[0 0; 0 9],'absolute',1e-10);
    assertEqual(dijInterval.OARRadiusRank,[1;1]);

function test_interval3_dynamic_radius_factor_respects_kmax
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    scenarioValues = zeros(3,5);
    pln.multScen = fixtureScenarioModel(ct,[1;1;1],scenarioValues, ...
        [1 1 1; 1 2 1; 1 3 1],true(1,3,1),ones(3,1)/3);
    dij.physicalDose = cell(1,3,1);
    dij.physicalDose{1} = sparse([1 0; 0 2; 1 0; 0 1]);
    dij.physicalDose{2} = sparse([1 0; 0 2; 0 1; 0 1]);
    dij.physicalDose{3} = sparse([1 0; 0 2; 2 2; 0 1]);
    cfg.KMode = 'dynamic';
    cfg.KMax = 1;
    cfg.RetentionThreshold = 1.0;

    [plnOut,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertEqual(dijInterval.OARRadiusRank(1),1);
    assertEqual(size(dijInterval.OARRadiusFactor{1},1),2);
    assertEqual(size(dijInterval.OARRadiusFactor{1},2),1);

function test_interval3_oar_radius_factor_accepts_sufficient_memory_limit
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.KMode = 'static';
    cfg.KMax = 2;
    [plnDefault,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijIntervalDefault = plnDefault.propOpt.dij_interval;

    cfg.MemoryLimitMB = 1;
    [plnLimited,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijIntervalLimited = plnLimited.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijIntervalLimited.center), ...
        full(dijIntervalDefault.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(reconstructOARRadiusCovariance(dijIntervalLimited,1)), ...
        full(reconstructOARRadiusCovariance(dijIntervalDefault,1)),'absolute',1e-12);

function test_interval3_oar_radius_factor_memory_scales_with_scenario_matrix
    [ct,cst,pln,~,cfg] = singleCtFixture();
    numBixels = 50;
    dij = baseDij(ct.cubeDim,numBixels);
    dij.physicalDose = cell(1,2,1);
    dij.physicalDose{1} = sparse(4,numBixels);
    dij.physicalDose{2} = sparse(4,numBixels);
    dij.physicalDose{1}(3,1:25) = 1;
    dij.physicalDose{2}(3,1:25) = 2;
    dij.physicalDose{1}(4,26:50) = 1;
    dij.physicalDose{2}(4,26:50) = 3;
    cfg.KMode = 'static';
    cfg.KMax = 2;
    cfg.MemoryLimitMB = 0.01;

    [plnOut,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    scenarioMatrix = [dij.physicalDose{1}(3,:); dij.physicalDose{2}(3,:)];
    centerRow = dijInterval.center(3,:);
    expectedCovariance = scenarioMatrix' * spdiags(dijInterval.scenarioWeights,0,2,2) * ...
        scenarioMatrix - centerRow' * centerRow;
    reconstructedCovariance = reconstructOARRadiusCovariance(dijInterval,1);

    assertEqual(size(dijInterval.OARRadiusFactor{1},1),numBixels);
    assertElementsAlmostEqual(full(reconstructedCovariance), ...
        full(expectedCovariance),'absolute',1e-10);

function test_interval3_dynamic_radius_factor_rank_selection_handles_large_energy_values
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.KMode = 'dynamic';
    cfg.KMax = 2;
    cfg.RetentionThreshold = 1.0;
    largeValue = 1e150;
    dij.physicalDose{1}(3,:) = sparse([largeValue 0]);
    dij.physicalDose{2}(3,:) = sparse([-largeValue 0]);
    dij.physicalDose{1}(4,:) = sparse([0 0]);
    dij.physicalDose{2}(4,:) = sparse([0 0]);

    [plnOut,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertEqual(dijInterval.OARRadiusRank(1),1);
    assertTrue(all(isfinite(nonzeros(dijInterval.OARRadiusFactor{1}))));

function test_interval3_oar_radius_factor_rejects_low_memory_limit
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.MemoryLimitMB = 1e-6;

    assertExceptionThrown(@() matRad_calcDoseInterval3(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_interval2_low_memory_limit_does_not_apply_oar_radius_factor_guard
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.MemoryLimitMB = 1e-6;

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijInterval.center(1:4,:)), ...
        [2.5 0; 0 3.5; 1.75 1; 0 4],'absolute',1e-12);

function test_interval3_batch_size_does_not_change_result_when_memory_allows
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.KMode = 'static';
    cfg.KMax = 2;
    cfg.MemoryLimitMB = 1;
    cfg.BatchSize = 1;
    [plnBatch,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijIntervalBatch = plnBatch.propOpt.dij_interval;

    cfg.BatchSize = 99;
    [plnFull,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijIntervalFull = plnFull.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijIntervalBatch.center), ...
        full(dijIntervalFull.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(reconstructOARRadiusCovariance(dijIntervalBatch,2)), ...
        full(reconstructOARRadiusCovariance(dijIntervalFull,2)),'absolute',1e-12);

function test_interval3_parallel_oar_center_radius_factor_matches_serial
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.KMode = 'static';
    cfg.KMax = 2;
    cfg.MemoryLimitMB = 1;
    cfg.BatchSize = 1;

    cfg.UseParallel = false;
    [plnSerial,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijIntervalSerial = plnSerial.propOpt.dij_interval;

    cfg.UseParallel = true;
    [plnParallel,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijIntervalParallel = plnParallel.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijIntervalParallel.center), ...
        full(dijIntervalSerial.center),'absolute',1e-12);
    assertElementsAlmostEqual(dijIntervalParallel.OARRadiusRank, ...
        dijIntervalSerial.OARRadiusRank, ...
        'absolute',1e-12);
    for i = 1:numel(dijIntervalSerial.OARRadiusFactor)
        serialCovariance = reconstructOARRadiusCovariance(dijIntervalSerial,i);
        parallelCovariance = reconstructOARRadiusCovariance(dijIntervalParallel,i);
        assertElementsAlmostEqual(full(parallelCovariance),full(serialCovariance), ...
            'absolute',1e-12);
    end

function test_interval3_collect_timing_reports_oar_radius_factor_components
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.KMode = 'static';
    cfg.KMax = 2;
    cfg.CollectTiming = true;
    cfg.UseParallel = false;

    [plnOut,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    timing = plnOut.propOpt.dij_interval.timing;

    assertEqual(timing.intervalMode,'INTERVAL3');
    assertFalse(timing.oar.parallelEnabled);
    assertTimingStageIsValid(timing.target);
    assertTimingStageIsValid(timing.oar);
    assertTrue(timing.target.radiusMultiplySeconds >= 0);
    assertTrue(timing.oar.factorSeconds >= 0);
    assertEqual(timing.oar.radiusMultiplySeconds,0);
    assertTrue(timing.oar.serialAssemblySeconds >= 0);

function test_interval3_zero_oar_covariance_is_valid_zero_rank
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    dij.physicalDose{2}(3,:) = dij.physicalDose{1}(3,:);
    cfg.OARStructSel = 'OAR';
    cfg.targetStructSel = 'PTV';

    [plnOut,~] = calcInterval3(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertEqual(dijInterval.OARRadiusRank(1),0);
    assertEqual(size(dijInterval.OARRadiusFactor{1}),[2 0]);

function test_interval2_multict_defaults_to_first_ct_scenario
    [ct,cst,pln,dij,cfg] = multiCtFixture(1);

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijInterval.center(:,1)), ...
        full(dij.physicalDose{1}),'absolute',1e-12);
    assertEqual(dijInterval.scenarioCtScenIds,1);
    assertElementsAlmostEqual(dijInterval.scenarioWeights,1,'absolute',1e-12);

    pln.propOpt.scen4D = [];
    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijInterval.center(:,1)), ...
        full(dij.physicalDose{1}),'absolute',1e-12);
    assertEqual(dijInterval.scenarioCtScenIds,1);
    assertElementsAlmostEqual(dijInterval.scenarioWeights,1,'absolute',1e-12);

function test_interval2_multict_all_maps_to_reference_scenario
    [ct,cst,pln,dij,cfg,expectedCenter] = multiCtFixture(1);
    pln.propOpt.scen4D = 'all';

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijInterval.center(:,1)),expectedCenter,'absolute',1e-12);
    assertEqual(dijInterval.refScen,1);
    assertEqual(dijInterval.scenarioCtScenIds,[1;2]);

function test_interval2_scen4d_filters_multict_scenarios
    [ct,cst,pln,dij,cfg] = multiCtFixture(1);
    pln.propOpt.scen4D = 1;

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijInterval.center(:,1)), ...
        full(dij.physicalDose{1}),'absolute',1e-12);
    assertEqual(dijInterval.scenarioCtScenIds,1);
    assertElementsAlmostEqual(dijInterval.scenarioWeights,1,'absolute',1e-12);

function test_interval2_scen4d_filters_and_maps_selected_ct_scenario
    [ct,cst,pln,dij,cfg] = multiCtFixture(1);
    pln.propOpt.scen4D = 2;

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    [xGrid,~,~] = meshgrid(1:ct.cubeDim(2),1:ct.cubeDim(1),1:ct.cubeDim(3));
    mappedRows = zeros(ct.cubeDim);
    mappedRows(:,2:end,:) = xGrid(:,1:end-1,:);
    assertElementsAlmostEqual(full(dijInterval.center(:,1)), ...
        mappedRows(:),'absolute',1e-12);
    assertEqual(dijInterval.scenarioCtScenIds,2);
    assertElementsAlmostEqual(dijInterval.scenarioWeights,1,'absolute',1e-12);

function test_interval2_scen4d_rejects_inactive_ct_scenario
    [ct,cst,pln,dij,cfg] = multiCtFixture(1);
    pln.propOpt.scen4D = 3;

    assertExceptionThrown(@() matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_interval2_multict_supports_nonfirst_reference_scenario
    [ct,cst,pln,dij,cfg,expectedCenter] = multiCtFixture(2);
    pln.propOpt.scen4D = 'all';

    [plnOut,dijIntervalContext] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijInterval.center(:,1)),expectedCenter,'absolute',1e-12);
    assertEqual(dijInterval.refScen,2);
    assertEqual(dijInterval.scenarioCtScenIds,[2;3]);
    assertEqual(dijIntervalContext.scenarioModel.numScenarios(),1);
    assertEqual(dijIntervalContext.scenarioModel.getDijScenarioIndex(1),1);
    assertEqual(dijIntervalContext.scenarioModel.getCtScenario(1),2);
    assertEqual(plnOut.multScen.getCtScenario(1),2);

function test_interval2_multict_requires_pull_dvf
    [ct,cst,pln,dij,cfg] = multiCtFixture(1);
    pln.propOpt.scen4D = 'all';
    ct.dvf = {};

    assertExceptionThrown(@() matRad_calcDoseInterval2(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_interval2_multict_uses_dij_ct_grid_when_ct_axes_are_missing
    [ct,cst,pln,dij,cfg,expectedCenter] = multiCtResamplingFixture();
    pln.propOpt.scen4D = 'all';

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);
    dijInterval = plnOut.propOpt.dij_interval;

    assertElementsAlmostEqual(full(dijInterval.center(:,1)),expectedCenter,'absolute',1e-12);

function test_interval2_multict_maps_y_displacement_and_boundaries
    [ct,cst,pln,dij,cfg,expectedCenter] = multiCtYDisplacementFixture();
    pln.propOpt.scen4D = 'all';

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(plnOut.propOpt.dij_interval.center(:,1)), ...
        expectedCenter,'absolute',1e-12);

function test_interval2_multict_converts_millimeter_dvf_to_dose_voxels
    [ct,cst,pln,dij,cfg,expectedCenter] = multiCtMillimeterDvfFixture();
    pln.propOpt.scen4D = 'all';

    [plnOut,~] = calcInterval2(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(plnOut.propOpt.dij_interval.center(:,1)), ...
        expectedCenter,'absolute',1e-12);

function test_interval2_streaming_inmemory_matches_existing_std
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    [plnLegacy,~] = calcInterval2(ct,cst,[],pln,dij,cfg);

    cfg.IntervalMode = 'INTERVAL2';
    [plnStream,dijStreamContext] = matRad_calcDoseIntervalStreaming(ct,cst,[],pln,dij,cfg);

    legacyInterval = plnLegacy.propOpt.dij_interval;
    streamInterval = plnStream.propOpt.dij_interval;
    assertElementsAlmostEqual(full(streamInterval.center), ...
        full(legacyInterval.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(streamInterval.radius), ...
        full(legacyInterval.radius),'absolute',1e-12);
    assertEqual(streamInterval.precomputeMode,'streaming');
    assertEqual(streamInterval.secondPassStrategy,'disk');
    assertEqual(dijStreamContext.numOfScenarios,1);

function test_interval2_streaming_partial_selection_matches_inmemory
    [ct,cst,pln,dij,cfg] = partialSelectionFixture();
    cfg.BatchSize = 1;
    [plnLegacy,~] = calcInterval2(ct,cst,[],pln,dij,cfg);

    [plnStreamDisk,~] = matRad_calcDoseInterval2Streaming(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(plnStreamDisk.propOpt.dij_interval.center), ...
        full(plnLegacy.propOpt.dij_interval.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(plnStreamDisk.propOpt.dij_interval.radius), ...
        full(plnLegacy.propOpt.dij_interval.radius),'absolute',1e-12);
    assertEqual(plnStreamDisk.propOpt.dij_interval.secondPassStrategy,'disk');

    cfg.SecondPassStrategy = 'recompute';
    [plnStreamRecompute,~] = matRad_calcDoseInterval2Streaming(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(plnStreamRecompute.propOpt.dij_interval.center), ...
        full(plnLegacy.propOpt.dij_interval.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(plnStreamRecompute.propOpt.dij_interval.radius), ...
        full(plnLegacy.propOpt.dij_interval.radius),'absolute',1e-12);
    assertEqual(plnStreamRecompute.propOpt.dij_interval.secondPassStrategy, ...
        'recompute');

function test_interval2_streaming_wrapper_does_not_require_interval_mode
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    [plnLegacy,~] = calcInterval2(ct,cst,[],pln,dij,cfg);

    [plnStream,~] = matRad_calcDoseInterval2Streaming(ct,cst,[],pln,dij,cfg);

    assertEqual(plnStream.propOpt.dij_interval.intervalMode,'INTERVAL2');
    assertElementsAlmostEqual(full(plnStream.propOpt.dij_interval.center), ...
        full(plnLegacy.propOpt.dij_interval.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(plnStream.propOpt.dij_interval.radius), ...
        full(plnLegacy.propOpt.dij_interval.radius),'absolute',1e-12);

function test_interval2_streaming_accepts_dij_without_cfg
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    [plnLegacy,~] = calcInterval2(ct,cst,[],pln,dij,cfg);

    [plnStream,~] = matRad_calcDoseInterval2Streaming(ct,cst,[],pln,dij);

    assertElementsAlmostEqual(full(plnStream.propOpt.dij_interval.center), ...
        full(plnLegacy.propOpt.dij_interval.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(plnStream.propOpt.dij_interval.radius), ...
        full(plnLegacy.propOpt.dij_interval.radius),'absolute',1e-12);

function test_interval2_streaming_recomputes_scenario_dij_without_dij_argument
    [ct,cst,pln,stf] = photonTestDataFixture();
    cfg.SecondPassStrategy = 'disk';
    cfg.KeepCache = false;
    cfg.BatchSize = 10000;
    cfg.targetStructSel = 1;

    [plnStream,dijStreamContext] = matRad_calcDoseInterval2Streaming(ct,cst,stf,pln,cfg);
    dijInterval = plnStream.propOpt.dij_interval;

    assertEqual(size(dijInterval.center,1),dijStreamContext.doseGrid.numOfVoxels);
    assertEqual(size(dijInterval.center,2),dijStreamContext.totalNumOfBixels);
    assertTrue(nnz(dijInterval.center) > 0);
    assertTrue(nnz(dijInterval.radius) > 0);
    assertFalse(isfield(dijInterval,'cacheDir'));
    assertStreamingSizeTotal(dijInterval);
    assertEqual(dijInterval.intervalMode,'INTERVAL2');
    assertEqual(dijStreamContext.numOfScenarios,1);

function test_interval_streaming_rejects_duplicate_precomputed_dij_inputs
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.IntervalMode = 'INTERVAL2';
    cfg.PrecomputedDij = dij;

    assertExceptionThrown(@() matRad_calcDoseIntervalStreaming(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_interval2_streaming_recompute_extreme_matches_existing
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.RadiusMode = 'extreme';
    [plnLegacy,~] = calcInterval2(ct,cst,[],pln,dij,cfg);

    cfg.IntervalMode = 'INTERVAL2';
    cfg.SecondPassStrategy = 'recompute';
    [plnStream,~] = matRad_calcDoseIntervalStreaming(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(plnStream.propOpt.dij_interval.center), ...
        full(plnLegacy.propOpt.dij_interval.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(plnStream.propOpt.dij_interval.radius), ...
        full(plnLegacy.propOpt.dij_interval.radius),'absolute',1e-12);
    assertStreamingSizeRecompute(plnStream.propOpt.dij_interval);

function test_interval3_streaming_wrapper_does_not_require_interval_mode
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.KMode = 'static';
    cfg.KMax = 2;
    [plnLegacy,~] = calcInterval3(ct,cst,[],pln,dij,cfg);

    [plnStream,~] = matRad_calcDoseInterval3Streaming(ct,cst,[],pln,dij,cfg);

    assertEqual(plnStream.propOpt.dij_interval.intervalMode,'INTERVAL3');
    assertElementsAlmostEqual(full(plnStream.propOpt.dij_interval.center), ...
        full(plnLegacy.propOpt.dij_interval.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(reconstructOARRadiusCovariance( ...
        plnStream.propOpt.dij_interval,1)), ...
        full(reconstructOARRadiusCovariance(plnLegacy.propOpt.dij_interval,1)), ...
        'absolute',1e-12);

function test_interval3_streaming_disk_cache_keeps_distinct_hash_folders
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.RadiusMode = 'extreme';
    cfg.KMode = 'static';
    cfg.KMax = 2;
    [plnLegacy,~] = calcInterval3(ct,cst,[],pln,dij,cfg);

    cacheRoot = tempname();
    cleanup = onCleanup(@() deleteDirIfExists(cacheRoot));
    cfg.IntervalMode = 'INTERVAL3';
    cfg.CacheRoot = cacheRoot;
    cfg.KeepCache = true;

    [plnStream1,~] = matRad_calcDoseIntervalStreaming(ct,cst,[],pln,dij,cfg);
    [plnStream2,~] = matRad_calcDoseIntervalStreaming(ct,cst,[],pln,dij,cfg);

    interval1 = plnStream1.propOpt.dij_interval;
    interval2 = plnStream2.propOpt.dij_interval;
    assertTrue(isfield(interval1,'cacheDir'));
    assertTrue(isfield(interval2,'cacheDir'));
    assertStreamingSizeDisk(interval1);
    assertStreamingSizeDisk(interval2);
    assertFalse(strcmp(interval1.cacheDir,interval2.cacheDir));
    assertEqual(exist(interval1.cacheDir,'dir'),7);
    assertEqual(exist(interval2.cacheDir,'dir'),7);
    assertEqual(exist(fullfile(interval1.cacheDir,'metadata.mat'),'file'),2);
    metadataData = load(fullfile(interval1.cacheDir,'metadata.mat'));
    assertTrue(isfield(metadataData,'metadata'));
    metadata = metadataData.metadata;
    assertEqual(metadata.calculationName,'dose interval');
    assertEqual(metadata.intervalMode,'INTERVAL3');
    assertEqual(metadata.radiusMode,'extreme');
    assertEqual(metadata.quantity,'physicalDose');
    assertEqual(metadata.refScen,1);
    assertEqual(metadata.scenarioDijIx,[1;2]);
    assertEqual(metadata.scenarioCtScenIds,[1;1]);
    assertElementsAlmostEqual(metadata.scenarioWeights,[0.25;0.75], ...
        'absolute',1e-12);
    assertTrue(~isempty(dir(fullfile(interval1.cacheDir,'scenario_*_oar_block_*.mat'))));

    legacyInterval = plnLegacy.propOpt.dij_interval;
    assertElementsAlmostEqual(full(interval1.center), ...
        full(legacyInterval.center),'absolute',1e-12);
    assertElementsAlmostEqual(full(interval1.radius), ...
        full(legacyInterval.radius),'absolute',1e-12);
    assertElementsAlmostEqual(full(reconstructOARRadiusCovariance(interval1,1)), ...
        full(reconstructOARRadiusCovariance(legacyInterval,1)),'absolute',1e-12);

function test_interval3_streaming_disk_cache_cleans_hash_folder_by_default
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cacheRoot = tempname();
    cleanup = onCleanup(@() deleteDirIfExists(cacheRoot));
    cfg.IntervalMode = 'INTERVAL3';
    cfg.CacheRoot = cacheRoot;
    cfg.KeepCache = false;
    cfg.KMode = 'static';
    cfg.KMax = 2;

    [plnStream,~] = matRad_calcDoseIntervalStreaming(ct,cst,[],pln,dij,cfg);

    dijInterval = plnStream.propOpt.dij_interval;
    assertFalse(isfield(dijInterval,'cacheDir'));
    assertStreamingSizeDisk(dijInterval);
    assertEqual(numel(listCacheRunDirs(cacheRoot)),0);

function test_interval2_streaming_multict_mapping_matches_existing
    [ct,cst,pln,dij,cfg,expectedCenter] = multiCtFixture(1);
    pln.propOpt.scen4D = 'all';
    cfg.IntervalMode = 'INTERVAL2';

    [plnStream,~] = matRad_calcDoseIntervalStreaming(ct,cst,[],pln,dij,cfg);

    assertElementsAlmostEqual(full(plnStream.propOpt.dij_interval.center(:,1)), ...
        expectedCenter,'absolute',1e-12);
    assertEqual(plnStream.propOpt.dij_interval.scenarioCtScenIds,[1;2]);

function test_interval_streaming_requires_interval_mode
    [ct,cst,pln,~,cfg] = singleCtFixture();

    assertExceptionThrown(@() matRad_calcDoseIntervalStreaming(ct,cst,[],pln,cfg), ...
        'matRad:Error');

function test_interval_streaming_rejects_invalid_second_pass_strategy
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.IntervalMode = 'INTERVAL2';
    cfg.SecondPassStrategy = 'memory';

    assertExceptionThrown(@() matRad_calcDoseIntervalStreaming(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function test_interval_streaming_rejects_invalid_cache_root
    [ct,cst,pln,~,cfg] = singleCtFixture();
    cfg.IntervalMode = 'INTERVAL2';
    cfg.CacheRoot = 1;

    assertExceptionThrown(@() matRad_calcDoseIntervalStreaming(ct,cst,[],pln,cfg), ...
        'matRad:Error');

function test_interval_streaming_rejects_cache_root_file
    [ct,cst,pln,dij,cfg] = singleCtFixture();
    cfg.IntervalMode = 'INTERVAL3';
    cfg.CacheRoot = tempname();
    fid = fopen(cfg.CacheRoot,'w');
    fwrite(fid,'not a folder');
    fclose(fid);
    cleanup = onCleanup(@() deleteFileIfExists(cfg.CacheRoot));

    assertExceptionThrown(@() matRad_calcDoseIntervalStreaming(ct,cst,[],pln,dij,cfg), ...
        'matRad:Error');

function covariance = reconstructOARRadiusCovariance(dijInterval,intervalIx)
    factor = dijInterval.OARRadiusFactor{intervalIx};
    covariance = factor * factor';

function [plnInterval,dijIntervalContext] = calcInterval2(ct,cst,stf,pln,dij,cfg)
    [plnInterval,dijIntervalContext] = matRad_calcDoseInterval2(ct,cst,stf,pln,dij,cfg);

function [plnInterval,dijIntervalContext] = calcInterval3(ct,cst,stf,pln,dij,cfg)
    [plnInterval,dijIntervalContext] = matRad_calcDoseInterval3(ct,cst,stf,pln,dij,cfg);

function assertTimingStageIsValid(stageTiming)
    assertTrue(isstruct(stageTiming));
    assertTrue(stageTiming.numVoxels >= 0);
    assertTrue(stageTiming.batchSize >= 0);
    assertTrue(stageTiming.numBatches >= 0);
    timingFields = {'extractMapSeconds','centerAccumSeconds', ...
        'radiusMultiplySeconds','factorSeconds','centeredRowsSeconds', ...
        'parallelSetupSeconds','parallelComputeWallSeconds', ...
        'serialAssemblySeconds','wallSeconds'};
    for fieldIx = 1:numel(timingFields)
        value = stageTiming.(timingFields{fieldIx});
        assertTrue(isnumeric(value));
        assertTrue(isscalar(value));
        assertTrue(isfinite(value));
        assertTrue(value >= 0);
    end

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
    scenMask = false(ct.numOfCtScen,1,1);
    scenMask(refScen) = true;
    scenMask(refScen + 1) = true;
    pln.multScen = fixtureScenarioModel(ct,[refScen; refScen + 1],zeros(2,5), ...
        [refScen 1 1; refScen + 1 1 1],scenMask,[0.5;0.5]);

    dij = baseDij(dim,1);
    dij.physicalDose = cell(ct.numOfCtScen,1,1);
    dij.physicalDose{refScen} = sparse(sourceRows);
    dij.physicalDose{refScen + 1} = sparse(sourceRows);
    cfg.refScen = refScen;
    cfg.BatchSize = 2;

    expectedCenter = 0.5*sourceRows + 0.5*mappedRows;

function [ct,cst,pln,dij,cfg,expectedCenter] = multiCtYDisplacementFixture()
    dim = [3 2 1];
    [~,yGrid,~] = meshgrid(1:dim(2),1:dim(1),1:dim(3));
    sourceRows = yGrid(:);
    mappedRows = zeros(dim);
    mappedRows(2:end,:,:) = yGrid(1:end-1,:,:);
    mappedRows = mappedRows(:);

    [ct,cst,pln,dij,cfg] = multiCtMappingBaseFixture(dim,sourceRows);
    ct.dvf{2}(2,:,:,:) = 1;

    expectedCenter = 0.5*sourceRows + 0.5*mappedRows;

function [ct,cst,pln,dij,cfg,expectedCenter] = multiCtMillimeterDvfFixture()
    dim = [2 3 1];
    [xGrid,~,~] = meshgrid(1:dim(2),1:dim(1),1:dim(3));
    sourceRows = xGrid(:);
    mappedRows = zeros(dim);
    mappedRows(:,2:end,:) = xGrid(:,1:end-1,:);
    mappedRows = mappedRows(:);

    [ct,cst,pln,dij,cfg] = multiCtMappingBaseFixture(dim,sourceRows);
    ct.resolution = struct('x',2,'y',1,'z',1);
    ct.dvfMetadata.dvfUnits = 'mm';
    ct.dvf{2}(1,:,:,:) = 2;
    dij.doseGrid.resolution = struct('x',2,'y',1,'z',1);
    dij.ctGrid = dij.doseGrid;

    expectedCenter = 0.5*sourceRows + 0.5*mappedRows;

function [ct,cst,pln,dij,cfg] = multiCtMappingBaseFixture(dim,sourceRows)
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

    cst = cell(1,6);
    cst = addStructure(cst,1,'PTV','TARGET',(1:prod(dim))');
    cst{1,4} = cell(1,2);
    cst{1,4}{1} = (1:prod(dim))';

    pln.bioParam.quantityOpt = 'physicalDose';
    pln.multScen = fixtureScenarioModel(ct,[1;2],zeros(2,5), ...
        [1 1 1; 2 1 1],true(2,1,1),[0.5;0.5]);

    dij = baseDij(dim,1);
    dij.physicalDose = cell(2,1,1);
    dij.physicalDose{1} = sparse(sourceRows);
    dij.physicalDose{2} = sparse(sourceRows);
    cfg.refScen = 1;
    cfg.BatchSize = 2;

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
    pln.multScen = fixtureScenarioModel(ct,[1;2],zeros(2,5), ...
        [1 1 1; 2 1 1],true(2,1,1),[0.5;0.5]);

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

function scenarioModel = fixtureScenarioModel(ct,ctScenIds,scenarioValues,linearMask,scenMask,scenarioWeights)
    scenarioModel = matRad_NominalScenario(ct);
    dimensions = matRad_createScenarioComponents([1 1 1],1,1);
    scenForProb = [ctScenIds(:) scenarioValues];
    scenarioModel.setScenarioRealizations(dimensions,scenarioValues,ctScenIds, ...
        scenarioWeights(:),scenarioWeights(:),scenForProb,linearMask,scenMask);

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
