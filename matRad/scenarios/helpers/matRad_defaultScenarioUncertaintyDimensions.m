function uncertaintyDimensionNames = matRad_defaultScenarioUncertaintyDimensions()
% matRad_defaultScenarioUncertaintyDimensions returns public uncertainty dimensions.
%
% call
%   uncertaintyDimensionNames = matRad_defaultScenarioUncertaintyDimensions()
%
% output
%   uncertaintyDimensionNames: ordered cell array of supported public
%                       uncertainty dimensions
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

dimensions = matRad_getScenarioDimensionRegistry([1 1 1], 1, 1, 0, 1, 1);
uncertaintyDimensionNames = {dimensions.name};

end
