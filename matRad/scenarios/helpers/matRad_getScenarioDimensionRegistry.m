function dimensions = matRad_getScenarioDimensionRegistry(shiftSD, rangeAbsSD, rangeRelSD, numOfBeams, gantryAngleSD, couchAngleSD)
% matRad_getScenarioDimensionRegistry returns registered uncertainty dimensions.
%
% call
%   dimensions = matRad_getScenarioDimensionRegistry(shiftSD,rangeAbsSD,rangeRelSD,numOfBeams,gantryAngleSD,couchAngleSD)
%
% input
%   shiftSD:             1x3 setup standard deviations [mm]
%   rangeAbsSD:          absolute range standard deviation [mm]
%   rangeRelSD:          relative range standard deviation [%]
%   numOfBeams:          number of beams used for beam-specific dimensions
%   gantryAngleSD:       gantry angle standard deviation [deg]
%   couchAngleSD:        couch angle standard deviation [deg]
%
% output
%   dimensions:          struct array describing uncertainty dimensions,
%                       components, units, scales, applicators, and context
%                       requirements
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

if nargin < 4 || isempty(numOfBeams)
    numOfBeams = 0;
end

if nargin < 5 || isempty(gantryAngleSD)
    gantryAngleSD = 1;
end

if nargin < 6 || isempty(couchAngleSD)
    couchAngleSD = 1;
end

dimensions = struct('name', {}, 'componentNames', {}, 'units', {}, ...
                    'defaultScale', {}, 'applicator', {}, 'requiresContext', {});

dimensions(end + 1) = helper_createDimension('ct', {}, {}, [], ...
                                             'matRad_CtScenarioApplicator', false);
dimensions(end + 1) = helper_createDimension('setup', ...
                                             {'setup.x', 'setup.y', 'setup.z'}, {'mm', 'mm', 'mm'}, shiftSD, ...
                                             'matRad_SetupShiftApplicator', false);
dimensions(end + 1) = helper_createDimension('range', ...
                                             {'range.absolute', 'range.relative'}, {'mm', 'fraction'}, ...
                                             [rangeAbsSD rangeRelSD ./ 100], 'matRad_RangeShiftApplicator', false);
dimensions(end + 1) = helper_createDimension('gantry', ...
                                             helper_beamComponentNames('gantry', numOfBeams), ...
                                             helper_repeatedUnits('deg', numOfBeams), ...
                                             gantryAngleSD .* ones(1, numOfBeams), ...
                                             'matRad_GantryAngleApplicator', true);
dimensions(end + 1) = helper_createDimension('couch', ...
                                             helper_beamComponentNames('couch', numOfBeams), ...
                                             helper_repeatedUnits('deg', numOfBeams), ...
                                             couchAngleSD .* ones(1, numOfBeams), ...
                                             'matRad_CouchAngleApplicator', true);

end

function dimension = helper_createDimension(name, componentNames, units, defaultScale, applicator, requiresContext)

dimension = struct();
dimension.name = name;
dimension.componentNames = componentNames;
dimension.units = units;
dimension.defaultScale = defaultScale;
dimension.applicator = applicator;
dimension.requiresContext = requiresContext;

end

function componentNames = helper_beamComponentNames(prefix, numOfBeams)

componentNames = cell(1, numOfBeams);
for i = 1:numOfBeams
    componentNames{i} = sprintf('%s.beam%d', prefix, i);
end

end

function units = helper_repeatedUnits(unit, numOfBeams)

units = repmat({unit}, 1, numOfBeams);

end
