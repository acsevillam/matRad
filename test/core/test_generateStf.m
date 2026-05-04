function test_suite = test_generateStf
    %The output should always be test_suite, and the function name the same as
    %your file name
    
    %% Header
    % The header is required to be in this format for automatic test collection
    % by MOxUnit
    
    %To collect all tests defined below, this is needed in newer Matlab
    %versions. test_functions will collect function handles to below test
    %functions
    test_functions=localfunctions(); 
    
    % This will initialize the test suite, i.e., take the functions from
    % test_functions, check if they contain "test", convert them into a MOxUnit
    % Test Case, and add them to the test-runner
    initTestSuite;
    
    %% Custom Tests
    % Tests use assert*-like Functions to check outputs etc:
    % assertTrue(a) - a is true
    % assertFalse(a) - a is false
    % assertEqual(a,b) - a and be are equal (isequal)
    % assertElementsAlmostEqual(a,b) - numerical test for all vector / matrix elements. Has Additional arguments for absolute / relative tolerance 
    % assertVectorsAlmostEqual(a,b) - numerical test using vector norm
    % assertExceptionThrown(f,id) - test if exception of id is thrown. Take care of Octave issues with exception id (or don't provide id)
    % Check MOxUnit for more information or look at other tests
    
    function test_generateStf_photons
        load TG119.mat;
        pln = helper_basicPln('photons',ct,cst);
        
        stf = matRad_generateStf(ct,cst,pln);
        assertTrue(isstruct(stf));
        assertEqual(numel(stf),numel(pln.propStf.gantryAngles));

    function test_generateStf_protons
        load TG119.mat;
        pln = helper_basicPln('protons',ct,cst);
        
        stf = matRad_generateStf(ct,cst,pln);
        assertTrue(isstruct(stf));
        assertEqual(numel(stf),numel(pln.propStf.gantryAngles));

    function test_generateStf_protonsAngularRandomScenarios
        load TG119.mat;
        pln = helper_basicPln('protons',ct,cst);
        pln.propStf.addMargin = false;
        pln.propOpt.runDAO = false;

        nominalStf = matRad_generateStf(ct,cst,pln);
        nominalRayPos = vertcat(nominalStf(1).ray.rayPos_bev);

        pln.multScen = matRad_RandomScenarios(ct);
        pln.multScen.nSamples = 2;
        pln.multScen.randomSeed = 13;
        pln.multScen.gantryAngleSD = 10;
        pln.multScen.couchAngleSD = 10;
        pln.multScen.scenarioDimensionActive = {'ct','gantry','couch'};

        stf = matRad_generateStf(ct,cst,pln);
        angularRayPos = vertcat(stf(1).ray.rayPos_bev);
        expectedRayPos = helper_projectTargetRayPositionsForScenarios(ct,cst,pln,1);

        assertTrue(isstruct(stf));
        assertEqual(numel(stf),numel(pln.propStf.gantryAngles));
        assertTrue(stf(1).numOfRays >= nominalStf(1).numOfRays);
        assertTrue(stf(1).totalNumOfBixels > 0);
        angularOnlyExpectedRayPos = expectedRayPos(~ismember(expectedRayPos,nominalRayPos,'rows'),:);
        assertFalse(isempty(angularOnlyExpectedRayPos));
        assertTrue(any(ismember(angularOnlyExpectedRayPos,angularRayPos,'rows')));

    function test_generateStf_helium
        load TG119.mat;
        pln = helper_basicPln('helium',ct,cst);
        
        stf = matRad_generateStf(ct,cst,pln);
        assertTrue(isstruct(stf));
        assertEqual(numel(stf),numel(pln.propStf.gantryAngles)); 

    function test_generateStf_carbon
        load TG119.mat;
        pln = helper_basicPln('carbon',ct,cst);
        
        stf = matRad_generateStf(ct,cst,pln);
        assertTrue(isstruct(stf));
        assertEqual(numel(stf),numel(pln.propStf.gantryAngles)); 

    function test_generateStf_noTargetObjectives
        load TG119.mat;
        [cst{:,6}] = deal([]);
        pln = helper_basicPln('photons',ct,cst);
                
        stf = matRad_generateStf(ct,cst,pln);
        assertTrue(isstruct(stf));
        assertEqual(numel(stf),numel(pln.propStf.gantryAngles));

    
    function pln = helper_basicPln(modality,ct,cst)
        pln.radiationMode = modality;
        pln.numOfFractions = 30;
        pln.machine = 'Generic';
        pln.multScen = matRad_NominalScenario(ct);
        pln.propStf.gantryAngles = 0;
        pln.propStf.couchAngles = 0;
        pln.propStf.bixelWidth = 5;
        pln.propStf.isoCenter = matRad_getIsoCenter(cst,ct);
        

    function rayPos = helper_projectTargetRayPositionsForScenarios(ct,cst,pln,beamIx)
        isTarget = cellfun(@(voiType) isequal(voiType,'TARGET'),cst(:,3));
        hasObjective = ~cellfun(@isempty,cst(:,6));
        useTargetForBixelPlacement = isTarget & hasObjective;
        if ~any(useTargetForBixelPlacement)
            useTargetForBixelPlacement(isTarget) = true;
        end

        V = [];
        for i = 1:size(cst,1)
            if useTargetForBixelPlacement(i)
                V = [V; cst{i,4}{1}];
            end
        end
        V = unique(V);

        ct = matRad_getWorldAxes(ct);
        voxTargetWorldCoords = matRad_cubeIndex2worldCoords(V,ct);
        isoCoords = voxTargetWorldCoords - pln.propStf.isoCenter(beamIx,:);
        scenarioIds = pln.multScen.scenarioIds();
        SAD = nominalSAD(pln);
        rayPos = zeros(0,3);

        for scenarioIx = 1:numel(scenarioIds)
            scenarioId = scenarioIds(scenarioIx);
            gantryOffsets = pln.multScen.getGantryAngleOffset(scenarioId);
            couchOffsets = pln.multScen.getCouchAngleOffset(scenarioId);
            gantryAngle = pln.propStf.gantryAngles(beamIx) + helper_beamAngleOffset(gantryOffsets,beamIx);
            couchAngle = pln.propStf.couchAngles(beamIx) + helper_beamAngleOffset(couchOffsets,beamIx);
            rotCoords = isoCoords * matRad_getRotationMatrix(gantryAngle,couchAngle);
            coordsAtIsoCenterPlane(:,1) = (rotCoords(:,1) * SAD) ./ (SAD + rotCoords(:,2));
            coordsAtIsoCenterPlane(:,2) = (rotCoords(:,3) * SAD) ./ (SAD + rotCoords(:,2));
            rayPos = [rayPos; unique(pln.propStf.bixelWidth * round([ ...
                coordsAtIsoCenterPlane(:,1) ...
                zeros(size(coordsAtIsoCenterPlane,1),1) ...
                coordsAtIsoCenterPlane(:,2)] / pln.propStf.bixelWidth),'rows')];
            coordsAtIsoCenterPlane = [];
        end

        rayPos = unique(rayPos,'rows');

    function SAD = nominalSAD(pln)
        machine = matRad_loadMachine(pln);
        SAD = machine.meta.SAD;

    function angleOffset = helper_beamAngleOffset(angleOffsets,beamIx)
        if isempty(angleOffsets) || beamIx > numel(angleOffsets)
            angleOffset = 0;
        else
            angleOffset = angleOffsets(beamIx);
        end
