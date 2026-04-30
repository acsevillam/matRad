function test_suite = test_compareDose

test_functions = localfunctions();

initTestSuite;

function test_compareDose_accepts_empty_cst_for_gamma_title
    oldFigureVisibility = get(0,'DefaultFigureVisible');
    cleanup = onCleanup(@() restoreFigureState(oldFigureVisibility));
    set(0,'DefaultFigureVisible','off');

    ct.cubeDim = [3 3 3];
    ct.numOfCtScen = 1;
    ct.resolution = struct('x',1,'y',1,'z',1);
    ct.cubeHU = {zeros(ct.cubeDim)};
    ct.x = 1:3;
    ct.y = 1:3;
    ct.z = 1:3;

    cube1 = ones(ct.cubeDim);
    cube2 = cube1;
    pln = struct('radiationMode','photons');

    [gammaCube,gammaPassRate] = matRad_compareDose(cube1,cube2,ct,[],[1 0 0],false,pln,[3 3],0,'global');

    assertEqual(size(gammaCube),ct.cubeDim);
    assertTrue(isnumeric(gammaPassRate));
    assertElementsAlmostEqual(gammaPassRate,100);

function restoreFigureState(oldFigureVisibility)
    close all;
    set(0,'DefaultFigureVisible',oldFigureVisibility);
