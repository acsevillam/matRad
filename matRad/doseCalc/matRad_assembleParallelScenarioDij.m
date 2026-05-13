function dij = matRad_assembleParallelScenarioDij(scenarioDijs,scenarioIds,scenarioModel)
% matRad_assembleParallelScenarioDij assembles scenario-local dijs
%
% call
%   dij = matRad_assembleParallelScenarioDij(scenarioDijs,scenarioIds,scenarioModel)
%
% input
%   scenarioDijs:    cell array with one single-scenario dij per scenario
%   scenarioIds:     public scenario ids from the original scenario model
%   scenarioModel:   original multi-scenario model
%
% output
%   dij:             robust multi-scenario dij with influence matrices at
%                    DIJ scenario indices and biological metadata at CT
%                    scenario indices
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

% Public wrapper kept for tests and compatibility; implementation lives in
% private to keep parallel-scenario assembly as internal infrastructure.
dij = matRad_assembleParallelScenarioDijCore(scenarioDijs,scenarioIds,scenarioModel);

end
