function [pln] = matRad_loadBeams(run_config,pln,ct,cst)

patient=run_config.description;
setup_type=run_config.plan_beams;

switch patient
    case 'prostate'
        switch setup_type
            case '5F'
                pln.numOfFractions         = 39;
                pln.propStf.gantryAngles   = [0:72:359];
                pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
                pln.propStf.bixelWidth     = 5;
            case '9F'
                pln.numOfFractions         = 39;
                pln.propStf.gantryAngles   = [0:40:359];
                pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
                pln.propStf.bixelWidth     = 5;
        end

    case 'breast'
        switch setup_type
            case '5F'
                pln.numOfFractions         = 16;
                pln.propStf.gantryAngles   = [0 45 90 135 315];
                pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
                pln.propStf.bixelWidth     = 5;
            case '7F'
                pln.numOfFractions         = 16;
                pln.propStf.gantryAngles   = [15 45 75 105 135 315 345];
                pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
                pln.propStf.bixelWidth     = 5;
        end

end

% Obtain the number of beams and voxels from the existing variables and
% calculate the iso-center which is per default the center of gravity of
% all target voxels.
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
    
end

