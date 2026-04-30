function test_suite = test_samplingAnalysis

test_functions = localfunctions();

initTestSuite;

function test_partial_sampling_skips_whole_ct_gamma
    [ct,cst,pln,caSampRes,mSampDose,resultGUInomScen] = samplingFixture();

    [~,doseStat] = matRad_samplingAnalysis(ct,cst,pln,caSampRes,mSampDose,resultGUInomScen);

    assertEqual(doseStat.gammaAnalysis.status,'skippedPartialSampling');
    assertEqual(doseStat.gammaAnalysis.gammaPassRate,[]);
    assertTrue(all(isnan(doseStat.gammaAnalysis.gammaCube(:))));
    assertElementsAlmostEqual(doseStat.sampleCoverageFraction,2/5);

function [ct,cst,pln,caSampRes,mSampDose,resultGUInomScen] = samplingFixture()
    ct = struct();
    ct.cubeDim = [5 1 1];
    ct.resolution = struct('x',1,'y',1,'z',1);
    ct.refScen = 1;

    cst = cell(1,6);
    cst{1,1} = 1;
    cst{1,2} = 'TARGET';
    cst{1,3} = 'TARGET';
    cst{1,4}{1} = [1 2];
    cst{1,5} = struct('Visible',true,'Priority',1);
    cst{1,6} = {};

    pln = struct();
    pln.subIx = [1;2];
    pln.bioParam = struct('quantityVis','physicalDose','model','none');
    pln.multScen = struct('totNumScen',2,'scenWeight',[1;1]);

    caSampRes = repmat(emptySampleResult(),1,2);
    caSampRes(1).dvh(1).volumePoints = [100 0];
    caSampRes(2).dvh(1).volumePoints = [100 50];
    caSampRes(1).qi(1).referenceDose = 2;
    caSampRes(2).qi(1).referenceDose = 2;

    mSampDose = single([1 2;3 4]);
    resultGUInomScen = struct('physicalDose',ones(ct.cubeDim));

function sampleResult = emptySampleResult()
    sampleResult = struct();
    sampleResult.dvh = struct('name','TARGET','doseGrid',[0 1],'volumePoints',[100 0]);
    sampleResult.qi = struct('name','TARGET','referenceDose',2);
