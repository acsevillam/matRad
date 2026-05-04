function test_suite = test_calcCubes

test_functions = localfunctions();

initTestSuite;

function test_calcCubesBedWithoutAlphaDoseUsesRBExD
dij.numOfBeams = 1;
dij.beamNum = 1;
dij.doseGrid.dimensions = [2 2 1];
dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);
dij.physicalDose{1} = sparse([0; 2; 0; 1]);
dij.ax{1} = 0.5*ones(dij.doseGrid.numOfVoxels,1);
dij.bx{1} = 0.05*ones(dij.doseGrid.numOfVoxels,1);
dij.RBE = 1.1;

resultGUI = matRad_calcCubes(1,dij);

expectedDose = full(dij.physicalDose{1}) * dij.RBE;
expectedBed = zeros(dij.doseGrid.numOfVoxels,1);
ixWeighted = full(dij.physicalDose{1}) > 0.01 * max(full(dij.physicalDose{1}));
expectedBed(ixWeighted) = (dij.ax{1}(ixWeighted) .* expectedDose(ixWeighted) + ...
    dij.bx{1}(ixWeighted) .* expectedDose(ixWeighted).^2) ./ dij.ax{1}(ixWeighted);

assertTrue(isfield(resultGUI,'RBExD'));
assertTrue(isfield(resultGUI,'BED'));
assertElementsAlmostEqual(resultGUI.BED(:),expectedBed,'absolute',1e-12);

function test_calcCubesDoseStdUsesQuadraticAccumulation
dij.numOfBeams = 2;
dij.beamNum = [1; 2];
dij.doseGrid.dimensions = [2 1 1];
dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);
dij.physicalDose{1} = sparse([1 2; 3 4]);
dij.physicalDose_std{1} = sparse([3 4; 5 12]);
dij.doseToWater{1} = sparse([2 3; 4 5]);
dij.doseToWater_std{1} = sparse([6 8; 7 24]);
dij.physicalDose_batchStd{1} = sparse([1 2; 2 3]);
w = [2; 3];

resultGUI = matRad_calcCubes(w,dij);

expectedPhysicalDose = full(dij.physicalDose{1}) * w;
expectedPhysicalDoseStd = sqrt(full(dij.physicalDose_std{1}).^2 * (w.^2));
expectedDoseToWaterStd = sqrt(full(dij.doseToWater_std{1}).^2 * (w.^2));
expectedBatchStd = sqrt(full(dij.physicalDose_batchStd{1}).^2 * (w.^2));
expectedBeam1Std = sqrt(full(dij.physicalDose_std{1}(:,1)).^2 * (w(1)^2));
expectedBeam2Std = sqrt(full(dij.physicalDose_std{1}(:,2)).^2 * (w(2)^2));

assertElementsAlmostEqual(resultGUI.physicalDose(:),expectedPhysicalDose,'absolute',1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_std(:),expectedPhysicalDoseStd,'absolute',1e-12);
assertElementsAlmostEqual(resultGUI.doseToWater_std(:),expectedDoseToWaterStd,'absolute',1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_batchStd(:),expectedBatchStd,'absolute',1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_std_beam1(:),expectedBeam1Std,'absolute',1e-12);
assertElementsAlmostEqual(resultGUI.physicalDose_std_beam2(:),expectedBeam2Std,'absolute',1e-12);
