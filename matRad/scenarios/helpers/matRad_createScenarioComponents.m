function components = matRad_createScenarioComponents(shiftSD,rangeAbsSD,rangeRelSD,activeDimensionNames,numOfBeams,gantryAngleSD,couchAngleSD)
% matRad_createScenarioComponents creates uncertainty realization components.
%
% Components stay present for compatibility. Their active flag is controlled
% by their parent uncertainty dimension/applicator in activeDimensionNames.
% Setup scales are stored in mm, absolute range scales in mm, relative range
% scales as fractions, and gantry/couch scales in degree.
%
% call
%   components = matRad_createScenarioComponents(shiftSD,rangeAbsSD,rangeRelSD)
%   components = matRad_createScenarioComponents(shiftSD,rangeAbsSD,rangeRelSD,activeDimensionNames,numOfBeams,gantryAngleSD,couchAngleSD)
%
% input
%   shiftSD:             1x3 setup standard deviations [mm]
%   rangeAbsSD:          absolute range standard deviation [mm]
%   rangeRelSD:          relative range standard deviation [%]
%   activeDimensionNames: active public uncertainty dimensions; supported
%                       values are ct, setup, range, gantry, and couch
%   numOfBeams:          number of beams used for per-beam angular
%                       components
%   gantryAngleSD:       gantry angle standard deviation [deg]
%   couchAngleSD:        couch angle standard deviation [deg]
%
% output
%   components:          struct array with name, applicator, unit, active,
%                       and scale fields
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    activeDimensionNames = matRad_defaultScenarioActiveDimensions();
end

if nargin < 5 || isempty(numOfBeams)
    numOfBeams = 0;
end

if nargin < 6 || isempty(gantryAngleSD)
    gantryAngleSD = 1;
end

if nargin < 7 || isempty(couchAngleSD)
    couchAngleSD = 1;
end

activeDimensionNames = matRad_normalizeScenarioDimensionActive(activeDimensionNames);
angleActive = any(strcmp(activeDimensionNames,'gantry')) || any(strcmp(activeDimensionNames,'couch'));

if angleActive
    validNumOfBeams = isnumeric(numOfBeams) && isscalar(numOfBeams) && ...
        isfinite(numOfBeams) && round(numOfBeams) == numOfBeams && numOfBeams >= 1;
    if ~validNumOfBeams
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Active gantry/couch uncertainty dimensions require numOfBeams to be a positive integer.');
    end
end

components = struct('name',{},'applicator',{},'unit',{},'active',{},'scale',{});
components(end+1) = createComponent('setup.x','setup','mm',shiftSD(1),activeDimensionNames);
components(end+1) = createComponent('setup.y','setup','mm',shiftSD(2),activeDimensionNames);
components(end+1) = createComponent('setup.z','setup','mm',shiftSD(3),activeDimensionNames);
components(end+1) = createComponent('range.absolute','range','mm',rangeAbsSD,activeDimensionNames);
components(end+1) = createComponent('range.relative','range','fraction',rangeRelSD./100,activeDimensionNames);

if angleActive
    for i = 1:numOfBeams
        components(end+1) = createComponent(sprintf('gantry.beam%d',i),'gantry','deg',gantryAngleSD,activeDimensionNames);
    end

    for i = 1:numOfBeams
        components(end+1) = createComponent(sprintf('couch.beam%d',i),'couch','deg',couchAngleSD,activeDimensionNames);
    end
end

matRad_validateScenarioComponents(components);

end

function component = createComponent(name,applicatorName,unit,scale,activeDimensionNames)

component = struct();
component.name = name;
component.applicator = applicatorName;
component.unit = unit;
component.active = any(strcmp(activeDimensionNames,applicatorName));
component.scale = scale;

end
