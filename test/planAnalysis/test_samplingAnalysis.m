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

function test_samplingReportScenarioParametersUseScenarioModelApi
    ct = struct('numOfCtScen',2);
    multScen = matRad_WorstCaseScenarios(ct);
    multScen.scenarioDimensionActive = {'ct','setup','range'};
    multScen.shiftSD = [1 2 3];
    multScen.rangeAbsSD = 4;
    multScen.rangeRelSD = 5;

    reportParameters = matRad_getSamplingReportScenarioParameters(multScen);

    assertEqual(reportParameters.numOfShiftScen,multScen.totNumShiftScen);
    assertEqual(reportParameters.numOfRangeShiftScen,multScen.totNumRangeScen);
    assertEqual(reportParameters.ctScen,'true');
    assertEqual(reportParameters.shiftScen,'true');
    assertEqual(reportParameters.rangeScen,'true');
    assertEqual(reportParameters.shiftCombType,multScen.combinations);

function test_calcStudyUsesPatientMatFilePathAndLocalTemplate
    calcStudyPath = which('matRad_calcStudy');
    calcStudySource = fileread(calcStudyPath);
    templatePath = fullfile(fileparts(calcStudyPath),'main_template.tex');

    assertEqual(nargin('matRad_calcStudy'),-2);
    assertTrue(exist(templatePath,'file') == 2);
    assertTrue(~isempty(strfind(calcStudySource,'load(matPatientPath')));
    assertTrue(~isempty(strfind(calcStudySource,'fileparts(mfilename(''fullpath''))')));
    assertTrue(isempty(strfind(calcStudySource,'exist(''matPatientPath'',''file'')')));
    assertTrue(isempty(strfind(calcStudySource,'tools'',''samplingAnalysis')));

function test_gaussianOrbitSamplesEigMethodPreservesCovarianceOrientation
    mu = [1; -2; 0.5];
    SIGMA = [3 1 0.4; 1 2 0.2; 0.4 0.2 1];
    radius = 1.7;
    nFrames = 12;

    samples = matRad_getGaussianOrbitSamples(mu,SIGMA,nFrames,radius,'Method','eig');

    centeredSamples = bsxfun(@minus,samples,mu);
    mahalanobisDistance = sum(centeredSamples .* (SIGMA \ centeredSamples),1);
    assertEqual(size(samples),[numel(mu) nFrames]);
    assertElementsAlmostEqual(mahalanobisDistance,radius^2*ones(1,nFrames),'absolute',1e-10);

function test_setupStudyTemplateUsesScenarioModelApi
    templatePath = which('matRad_setupStudyTemplate');
    templateSource = fileread(templatePath);
    forbiddenTerms = {'numOfShiftScen','numOfRangeShiftScen','shiftGenType', ...
        'scenCombType','matRad_createValidInstance','setup.x','range.absolute'};

    for i = 1:numel(forbiddenTerms)
        assertTrue(isempty(strfind(templateSource,forbiddenTerms{i})));
    end

    assertTrue(~isempty(strfind(templateSource,'scenarioDimensionActive = {''ct'',''setup'',''range''}')));
    assertTrue(~isempty(strfind(templateSource,'matRad_calcStudy(multScen')));

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
