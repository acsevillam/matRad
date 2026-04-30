function test_suite = test_elasticImageRegistration

test_functions = localfunctions();

initTestSuite;

function test_calc_dvf_fails_when_image_processing_toolbox_is_missing
    if hasElasticRegistrationDependencies()
        moxunit_throw_test_skipped_exception('Image Processing Toolbox is available.');
    end

    [ct,cst] = registrationFixture();
    obj = matRad_ElasticImageRegistration(ct,cst,1,struct('dvfType','pull'));

    assertExceptionThrown(@() obj.calcDVF(),'matRad:Error');

function test_prop_contours_preserves_empty_structures
    if ~hasElasticRegistrationDependencies()
        moxunit_throw_test_skipped_exception('Image Processing Toolbox is not available.');
    end

    [ct,cst] = registrationFixture();
    ct.dvfType = 'push';
    ct.dvfUnits = 'voxel';
    ct.dvf = {zeros([3 ct.cubeDim]),zeros([3 ct.cubeDim])};

    obj = matRad_ElasticImageRegistration(ct,cst,1,struct('dvfType','push'));
    [~,cstOut] = obj.propContours();

    assertEqual(size(cstOut,1),1);
    assertTrue(isempty(cstOut{1,4}{1,1}));
    assertTrue(isempty(cstOut{1,4}{1,2}));

function [ct,cst] = registrationFixture()
    ct = struct();
    ct.numOfCtScen = 2;
    ct.cubeDim = [2 2 1];
    ct.cubeHU = {zeros(ct.cubeDim),zeros(ct.cubeDim)};
    ct.resolution = struct('x',1,'y',1,'z',1);

    cst = cell(1,6);
    cst{1,1} = 1;
    cst{1,2} = 'EMPTY';
    cst{1,3} = 'TARGET';
    cst{1,4}{1,1} = [];
    cst{1,5} = struct('Visible',true);
    cst{1,6} = {};

function tf = hasElasticRegistrationDependencies()
    env = matRad_getEnvironment();
    tf = strcmp(env,'MATLAB') && exist('imregdemons','file') == 2 && exist('imwarp','file') == 2;
    if tf && (exist('license','builtin') == 5 || exist('license','file') == 2)
        try
            tf = license('test','image_toolbox');
        catch
        end
    end
