function test_suite = test_mapDoseToReferenceScenario

test_functions = localfunctions();

initTestSuite;

function test_referenceScenarioReturnsInputDose
ct = helper_createCtWithPullDvf();
doseCube = helper_createColumnDoseCube(ct.cubeDim);

[doseCubeRef, meta] = matRad_mapDoseToReferenceScenario(ct, doseCube, 1, 1);

assertElementsAlmostEqual(doseCubeRef, doseCube, 'absolute', 1e-12);
assertFalse(meta.mapped);
assertEqual(meta.sourceCtScenId, 1);
assertEqual(meta.referenceCtScenId, 1);

function test_pullDvfMapsDoseToReferenceScenario
ct = helper_createCtWithPullDvf();
doseCube = helper_createColumnDoseCube(ct.cubeDim);

doseCubeRef = matRad_mapDoseToReferenceScenario(ct, doseCube, 2, 1);

expectedCube = zeros(ct.cubeDim);
expectedCube(:, 2:end, :) = doseCube(:, 1:end - 1, :);
assertElementsAlmostEqual(doseCubeRef, expectedCube, 'absolute', 1e-12);

function test_missingDvfThrowsError
ct = helper_createCtWithPullDvf();
ct.dvf{2} = [];
doseCube = helper_createColumnDoseCube(ct.cubeDim);

assertExceptionThrown(@() matRad_mapDoseToReferenceScenario(ct, doseCube, 2, 1), ...
                      'matRad:Error');

function test_nonPullDvfThrowsError
ct = helper_createCtWithPullDvf();
ct.dvfMetadata.dvfType = 'push';
doseCube = helper_createColumnDoseCube(ct.cubeDim);

assertExceptionThrown(@() matRad_mapDoseToReferenceScenario(ct, doseCube, 2, 1), ...
                      'matRad:Error');

function test_referenceMetadataMismatchThrowsError
ct = helper_createCtWithPullDvf();
ct.dvfMetadata.refScen = 2;
ct.dvfMetadata.referenceCtScen = 2;
doseCube = helper_createColumnDoseCube(ct.cubeDim);

assertExceptionThrown(@() matRad_mapDoseToReferenceScenario(ct, doseCube, 2, 1), ...
                      'matRad:Error');

function test_invalidDoseCubeSizeThrowsError
ct = helper_createCtWithPullDvf();
doseCube = helper_createColumnDoseCube(ct.cubeDim);

assertExceptionThrown(@() matRad_mapDoseToReferenceScenario(ct, doseCube(:, :, 1), 2, 1), ...
                      'matRad:Error');

function ct = helper_createCtWithPullDvf()
ct.cubeDim = [3 4 2];
ct.numOfCtScen = 2;
ct.refScen = 1;
ct.resolution.x = 1;
ct.resolution.y = 1;
ct.resolution.z = 1;
ct.dvfMetadata.dvfType = 'pull';
ct.dvfMetadata.dvfUnits = 'voxel';
ct.dvfMetadata.refScen = 1;
ct.dvfMetadata.referenceCtScen = 1;

ct.dvf = cell(1, 2);
ct.dvf{1} = zeros([3 ct.cubeDim]);
ct.dvf{2} = zeros([3 ct.cubeDim]);
ct.dvf{2}(1, :, :, :) = 1;

function doseCube = helper_createColumnDoseCube(cubeDim)
doseCube = zeros(cubeDim);
for x = 1:cubeDim(2)
    doseCube(:, x, :) = x;
end
