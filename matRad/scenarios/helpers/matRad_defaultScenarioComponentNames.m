function componentNames = matRad_defaultScenarioComponentNames(numOfBeams)
% matRad_defaultScenarioComponentNames returns realization component names.
%
% call
%   componentNames = matRad_defaultScenarioComponentNames()
%   componentNames = matRad_defaultScenarioComponentNames(numOfBeams)
%
% input
%   numOfBeams:          number of beams for per-beam gantry and couch
%                       components
%
% output
%   componentNames:      ordered component names used by scenarioValues
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

if nargin < 1 || isempty(numOfBeams)
    numOfBeams = 0;
end

componentNames = {'setup.x','setup.y','setup.z','range.absolute','range.relative'};

for i = 1:numOfBeams
    componentNames{end+1} = sprintf('gantry.beam%d',i);
end

for i = 1:numOfBeams
    componentNames{end+1} = sprintf('couch.beam%d',i);
end

end
