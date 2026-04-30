function test_suite = test_mapDoseToReferenceScenario

test_functions = localfunctions();
initTestSuite;

function test_mapDoseToReferenceScenarioLeavesReferenceDoseUnchanged()
    ct = createTestCt('voxel');
    doseCube = reshape(1:prod(ct.cubeDim),ct.cubeDim);

    doseCubeRef = matRad_mapDoseToReferenceScenario(ct,doseCube,1,1);

    assertEqual(doseCubeRef,doseCube);

function test_mapDoseToReferenceScenarioUsesPullDvfInVoxelUnits()
    ct = createTestCt('voxel');
    [xGrid,~,~] = meshgrid(1:ct.cubeDim(2),1:ct.cubeDim(1),1:ct.cubeDim(3));
    doseCube = xGrid;
    ct.dvf{2}(1,:,:,:) = 1;

    doseCubeRef = matRad_mapDoseToReferenceScenario(ct,doseCube,2,1);

    expected = zeros(ct.cubeDim);
    expected(:,2:end,:) = xGrid(:,1:end-1,:);
    assertElementsAlmostEqual(doseCubeRef,expected);

function test_mapDoseToReferenceScenarioUsesPullDvfInMillimeters()
    ct = createTestCt('mm');
    [xGrid,~,~] = meshgrid(1:ct.cubeDim(2),1:ct.cubeDim(1),1:ct.cubeDim(3));
    doseCube = xGrid;
    ct.dvf{2}(1,:,:,:) = ct.resolution.x;

    doseCubeRef = matRad_mapDoseToReferenceScenario(ct,doseCube,2,1);

    expected = zeros(ct.cubeDim);
    expected(:,2:end,:) = xGrid(:,1:end-1,:);
    assertElementsAlmostEqual(doseCubeRef,expected);

function test_mapDoseToReferenceScenarioSupportsNonFirstReferenceScenario()
    ct = createTestCt('voxel');
    ct.numOfCtScen = 3;
    ct.refScen = 2;
    ct.dvfMetadata.refScen = 2;
    ct.dvfMetadata.referenceCtScen = 2;
    ct.dvf{3} = zeros([3 ct.cubeDim]);
    [xGrid,~,~] = meshgrid(1:ct.cubeDim(2),1:ct.cubeDim(1),1:ct.cubeDim(3));
    doseCube = xGrid;
    ct.dvf{3}(1,:,:,:) = 1;

    [doseCubeRef,meta] = matRad_mapDoseToReferenceScenario(ct,doseCube,3,2);

    expected = zeros(ct.cubeDim);
    expected(:,2:end,:) = xGrid(:,1:end-1,:);
    assertElementsAlmostEqual(doseCubeRef,expected);
    assertEqual(meta.sourceCtScen,3);
    assertEqual(meta.referenceCtScen,2);
    assertTrue(meta.mapped);

function test_mapDoseToReferenceScenarioRejectsMismatchedReferenceScenario()
    ct = createTestCt('voxel');
    ct.refScen = 2;
    ct.dvfMetadata.refScen = 2;
    ct.dvfMetadata.referenceCtScen = 2;
    doseCube = zeros(ct.cubeDim);

    assertExceptionThrown(@() matRad_mapDoseToReferenceScenario(ct,doseCube,1,1),'matRad:Error');

function ct = createTestCt(dvfUnits)
    ct = struct();
    ct.numOfCtScen = 2;
    ct.cubeDim = [4 5 3];
    ct.resolution = struct('x',2,'y',3,'z',4);
    ct.refScen = 1;
    ct.dvfMetadata.dvfType = 'pull';
    ct.dvfMetadata.dvfUnits = dvfUnits;
    ct.dvfMetadata.refScen = 1;
    ct.dvfMetadata.referenceCtScen = 1;
    ct.dvf = cell(1,2);
    ct.dvf{1} = zeros([3 ct.cubeDim]);
    ct.dvf{2} = zeros([3 ct.cubeDim]);
