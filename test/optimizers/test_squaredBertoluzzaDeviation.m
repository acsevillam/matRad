function test_suite = test_squaredBertoluzzaDeviation

test_functions=localfunctions();

initTestSuite;

function test_squared_bertoluzza_construct_and_reject_nominal_dose_call

    obj = DoseObjectives.matRad_SquaredBertoluzzaDeviation(2,5);

    assertEqual(obj.name,'Squared Bertoluzza Deviation');
    assertEqual(obj.penalty,2);
    assertEqual(obj.parameters{1},5);

    assertExceptionThrown(@() obj.computeDoseObjectiveFunction([3;7]), ...
        'matRad:Error');
    assertExceptionThrown(@() obj.computeDoseObjectiveGradient([3;7]), ...
        'matRad:Error');

function test_squared_bertoluzza_interval_value_and_gradient

    obj = DoseObjectives.matRad_SquaredBertoluzzaDeviation(3,1);
    [w,Ix,theta,dij_interval] = buildIntervalInputs();

    f = obj.computeDoseObjectiveFunction(w,Ix,theta,dij_interval);
    g = obj.computeFluenceObjectiveGradient(w,Ix,theta,dij_interval);

    assertElementsAlmostEqual(f,8.1,'absolute',1e-12);
    assertElementsAlmostEqual(g,[2.4;11.4],'absolute',1e-10);

function test_squared_bertoluzza_gradient_matches_central_difference

    obj = DoseObjectives.matRad_SquaredBertoluzzaDeviation(3,1);
    [w,Ix,theta,dij_interval] = buildIntervalInputs();

    f = @(x) obj.computeDoseObjectiveFunction(x,Ix,theta,dij_interval);
    g = obj.computeFluenceObjectiveGradient(w,Ix,theta,dij_interval);
    gNum = centralDifference(f,w);

    assertElementsAlmostEqual(g,gNum,'absolute',1e-6);

function test_squared_bertoluzza_struct_initialization_and_robustness

    s.className = 'DoseObjectives.matRad_SquaredBertoluzzaDeviation';
    s.parameters = {42};
    s.penalty = 7;
    s.robustness = 'INTERVAL2';

    obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(s);

    assertTrue(isa(obj,'DoseObjectives.matRad_SquaredBertoluzzaDeviation'));
    assertEqual(obj.parameters{1},42);
    assertEqual(obj.penalty,7);
    assertEqual(obj.robustness,'INTERVAL2');
    assertEqual(obj.availableRobustness(),{'INTERVAL2','INTERVAL3'});
    assertExceptionThrown(@() setObjectiveRobustness(obj,'none'),'matRad:Error');

function test_squared_bertoluzza_interval_objective_wrapper

    [optiProb,dij,cst,w] = buildIntervalOptimizationProblem();

    f = optiProb.matRad_objectiveFunction(w,dij,cst);
    g = optiProb.matRad_objectiveGradient(w,dij,cst);

    assertElementsAlmostEqual(f,8.1,'absolute',1e-12);
    assertElementsAlmostEqual(g,[2.4;11.4],'absolute',1e-10);

function test_interval3_target_uses_bertoluzza_target_model

    [optiProb,dij,cst,w] = buildIntervalOptimizationProblem();
    cst{1,6}{1}.robustness = 'INTERVAL3';

    f = optiProb.matRad_objectiveFunction(w,dij,cst);
    g = optiProb.matRad_objectiveGradient(w,dij,cst);

    assertElementsAlmostEqual(f,8.1,'absolute',1e-12);
    assertElementsAlmostEqual(g,[2.4;11.4],'absolute',1e-10);

function test_interval2_target_uses_reference_scenario_contour

    [optiProb,dij,cst,w] = buildReferenceContourTargetProblem('INTERVAL2');

    f = optiProb.matRad_objectiveFunction(w,dij,cst);
    g = optiProb.matRad_objectiveGradient(w,dij,cst);

    assertElementsAlmostEqual(f,8.1,'absolute',1e-12);
    assertElementsAlmostEqual(g,[2.4;11.4],'absolute',1e-10);

function test_interval3_target_uses_reference_scenario_contour

    [optiProb,dij,cst,w] = buildReferenceContourTargetProblem('INTERVAL3');

    f = optiProb.matRad_objectiveFunction(w,dij,cst);
    g = optiProb.matRad_objectiveGradient(w,dij,cst);

    assertElementsAlmostEqual(f,8.1,'absolute',1e-12);
    assertElementsAlmostEqual(g,[2.4;11.4],'absolute',1e-10);

function test_interval2_oar_uses_center_dose

    [optiProb,dij,cst,w] = buildIntervalOARCenterProblem('INTERVAL2');

    f = optiProb.matRad_objectiveFunction(w,dij,cst);
    g = optiProb.matRad_objectiveGradient(w,dij,cst);

    assertElementsAlmostEqual(f,9,'absolute',1e-12);
    assertElementsAlmostEqual(g,[0;12],'absolute',1e-12);

function test_interval2_oar_uses_reference_scenario_contour

    [optiProb,dij,cst,w] = buildReferenceContourOARCenterProblem('INTERVAL2');

    f = optiProb.matRad_objectiveFunction(w,dij,cst);
    g = optiProb.matRad_objectiveGradient(w,dij,cst);

    assertElementsAlmostEqual(f,9,'absolute',1e-12);
    assertElementsAlmostEqual(g,[0;12],'absolute',1e-12);

function test_interval3_oar_uses_center_plus_radius

    [optiProb,dij,cst,w,expectedDose] = buildIntervalOARRadiusProblem();
    objective = cst{1,6}{1};

    f = optiProb.matRad_objectiveFunction(w,dij,cst);
    expectedF = objective.penalty*objective.computeDoseObjectiveFunction(expectedDose);

    assertElementsAlmostEqual(f,expectedF,'absolute',1e-12);

function test_interval3_oar_gradient_matches_central_difference

    [optiProb,dij,cst,w] = buildIntervalOARRadiusProblem();

    f = @(x) optiProb.matRad_objectiveFunction(x,dij,cst);
    g = optiProb.matRad_objectiveGradient(w,dij,cst);
    gNum = centralDifference(f,w);

    assertElementsAlmostEqual(g,gNum,'absolute',1e-6);

function test_interval3_oar_uses_reference_scenario_contour

    [optiProb,dij,cst,w,expectedDose] = buildReferenceContourOARRadiusProblem();
    objective = cst{1,6}{1};

    f = optiProb.matRad_objectiveFunction(w,dij,cst);
    g = optiProb.matRad_objectiveGradient(w,dij,cst);
    gNum = centralDifference(@(x) optiProb.matRad_objectiveFunction(x,dij,cst),w);
    expectedF = objective.penalty*objective.computeDoseObjectiveFunction(expectedDose);

    assertElementsAlmostEqual(f,expectedF,'absolute',1e-12);
    assertElementsAlmostEqual(g,gNum,'absolute',1e-6);

function test_general_objectives_accept_interval_options_and_gui_lists_them

    objective = DoseObjectives.matRad_SquaredDeviation(1,0);
    objective.robustness = 'INTERVAL2';
    assertEqual(objective.robustness,'INTERVAL2');

    objective.robustness = 'INTERVAL3';
    assertEqual(objective.robustness,'INTERVAL3');

    widgetText = fileread(which('matRad_OptimizationWidget'));
    assertTrue(~isempty(strfind(widgetText,'INTERVAL2')));
    assertTrue(~isempty(strfind(widgetText,'INTERVAL3')));

function test_interval_target_rejects_non_bertoluzza_objective

    [optiProb,dij,cst,w] = buildIntervalOptimizationProblem();
    objective = DoseObjectives.matRad_SquaredDeviation(1,1);
    objective.robustness = 'INTERVAL2';
    cst{1,6} = {objective};

    assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w,dij,cst),'matRad:Error');

function test_interval3_oar_rejects_missing_svd_fields

    [optiProb,dij,cst,w] = buildIntervalOARCenterProblem('INTERVAL3');

    assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w,dij,cst),'matRad:Error');

function test_interval3_oar_rejects_uncovered_voxel

    [optiProb,dij,cst,w] = buildIntervalOARRadiusProblem();
    optiProb.dij_interval.OARSubIx = 1;
    optiProb.dij_interval.U = optiProb.dij_interval.U(1);
    optiProb.dij_interval.S = optiProb.dij_interval.S(1);
    optiProb.dij_interval.V = optiProb.dij_interval.V(1);

    assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w,dij,cst),'matRad:Error');

function test_interval_rejects_quantity_mismatch

    [optiProb,dij,cst,w] = buildIntervalOptimizationProblem();
    optiProb.quantityOpt = 'RBExD';
    optiProb.dij_interval.quantity = 'physicalDose';

    assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w,dij,cst),'matRad:Error');

function test_interval3_oar_accepts_zero_rank_radius

    [optiProb,dij,cst,w] = buildIntervalOARRadiusProblem();
    optiProb.dij_interval.U{1} = sparse(2,0);
    optiProb.dij_interval.S{1} = sparse(0,0);
    optiProb.dij_interval.V{1} = sparse(2,0);

    f = optiProb.matRad_objectiveFunction(w,dij,cst);

    assertTrue(isfinite(f));

function test_interval_rejects_missing_reference_contour

    [optiProb,dij,cst,w] = buildIntervalOptimizationProblem();
    optiProb.dij_interval.refScen = 2;
    cst{1,4} = {[1;2]};

    assertExceptionThrown(@() optiProb.matRad_objectiveFunction(w,dij,cst), ...
        'matRad:Error');

function [w,Ix,theta,dij_interval] = buildIntervalInputs()

    w = [1;2];
    Ix = [1;2];
    theta = 0.4;
    dij_interval.center = sparse([1 0; 0 2; 3 0]);
    dij_interval.radius = sparse([2 0.5; 0.5 1]);

function [optiProb,dij,cst,w] = buildIntervalOptimizationProblem()

    [w,Ix,theta,dij_interval] = buildIntervalInputs();

    backProjection = matRad_DoseProjection();
    backProjection.scenarios = 1;
    backProjection.scenarioProb = 1;
    backProjection.nominalCtScenarios = 1;

    optiProb = matRad_OptimizationProblem(backProjection);
    optiProb.theta1 = theta;
    optiProb.dij_interval = dij_interval;

    dij.physicalDose = {sparse(3,2)};
    dij.doseGrid.numOfVoxels = 3;
    dij.totalNumOfBixels = 2;

    objective = DoseObjectives.matRad_SquaredBertoluzzaDeviation(3,1);
    objective.robustness = 'INTERVAL2';

    cst = cell(1,6);
    cst{1,3} = 'TARGET';
    cst{1,4} = {Ix};
    cst{1,5}.alphaX = [];
    cst{1,5}.betaX = [];
    cst{1,6} = {objective};

function [optiProb,dij,cst,w] = buildReferenceContourTargetProblem(robustness)

    [optiProb,dij,cst,w] = buildIntervalOptimizationProblem();
    optiProb.dij_interval.refScen = 2;
    optiProb.dij_interval.center = sparse([99 0; 1 0; 0 2]);
    cst{1,4} = {1; [2;3]};
    cst{1,6}{1}.robustness = robustness;

    backProjection = matRad_DoseProjection();
    backProjection.scenarios = [1;2];
    backProjection.scenarioProb = [0.5;0.5];
    backProjection.nominalCtScenarios = 1;
    optiProb.BP = backProjection;

    dij.physicalDose = cell(2,1,1);
    dij.physicalDose{1} = sparse(3,2);
    dij.physicalDose{2} = sparse(3,2);

function [optiProb,dij,cst,w] = buildIntervalOARCenterProblem(robustness)

    [optiProb,dij,cst,w] = buildIntervalOptimizationProblem();

    objective = DoseObjectives.matRad_SquaredDeviation(2,1);
    objective.robustness = robustness;

    cst{1,3} = 'OAR';
    cst{1,6} = {objective};

function [optiProb,dij,cst,w] = buildReferenceContourOARCenterProblem(robustness)

    [optiProb,dij,cst,w] = buildIntervalOARCenterProblem(robustness);
    optiProb.dij_interval.refScen = 2;
    optiProb.dij_interval.center = sparse([99 0; 1 0; 0 2]);
    cst{1,4} = {1; [2;3]};

    backProjection = matRad_DoseProjection();
    backProjection.scenarios = [1;2];
    backProjection.scenarioProb = [0.5;0.5];
    backProjection.nominalCtScenarios = 1;
    optiProb.BP = backProjection;

    dij.physicalDose = cell(2,1,1);
    dij.physicalDose{1} = sparse(3,2);
    dij.physicalDose{2} = sparse(3,2);

function [optiProb,dij,cst,w,expectedDose] = buildIntervalOARRadiusProblem()

    w = [3;4];

    backProjection = matRad_DoseProjection();
    backProjection.scenarios = 1;
    backProjection.scenarioProb = 1;
    backProjection.nominalCtScenarios = 1;

    optiProb = matRad_OptimizationProblem(backProjection);
    optiProb.theta2 = 0.5;
    optiProb.dij_interval.center = sparse([1 0; 0 1; 0 0]);
    optiProb.dij_interval.radius = sparse(2,2);
    optiProb.dij_interval.OARSubIx = [1;2];
    optiProb.dij_interval.U = {speye(2),speye(2)};
    optiProb.dij_interval.S = {sparse(diag([4 1])),sparse(diag([1 9]))};
    optiProb.dij_interval.V = {speye(2),speye(2)};

    dCenter = optiProb.dij_interval.center([1;2],:)*w;
    dRadius = [sqrt(w'*diag([4 1])*w); sqrt(w'*diag([1 9])*w)];
    expectedDose = dCenter + optiProb.theta2*dRadius;

    dij.physicalDose = {sparse(3,2)};
    dij.doseGrid.numOfVoxels = 3;
    dij.totalNumOfBixels = 2;

    objective = DoseObjectives.matRad_SquaredDeviation(2,1);
    objective.robustness = 'INTERVAL3';

    cst = cell(1,6);
    cst{1,3} = 'OAR';
    cst{1,4} = {[1;2]};
    cst{1,5}.alphaX = [];
    cst{1,5}.betaX = [];
    cst{1,6} = {objective};

function [optiProb,dij,cst,w,expectedDose] = buildReferenceContourOARRadiusProblem()

    [optiProb,dij,cst,w,expectedDose] = buildIntervalOARRadiusProblem();
    optiProb.dij_interval.refScen = 2;
    cst{1,4} = {3; [1;2]};

    backProjection = matRad_DoseProjection();
    backProjection.scenarios = [1;2];
    backProjection.scenarioProb = [0.5;0.5];
    backProjection.nominalCtScenarios = 1;
    optiProb.BP = backProjection;

    dij.physicalDose = cell(2,1,1);
    dij.physicalDose{1} = sparse(3,2);
    dij.physicalDose{2} = sparse(3,2);

function g = centralDifference(f,w)

    delta = 1e-6;
    g = zeros(size(w));

    for i = 1:numel(w)
        step = zeros(size(w));
        step(i) = delta;
        g(i) = (f(w + step) - f(w - step)) / (2*delta);
    end

function setObjectiveRobustness(obj,robustness)

    obj.robustness = robustness;
