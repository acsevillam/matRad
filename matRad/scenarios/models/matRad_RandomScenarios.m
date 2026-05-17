classdef matRad_RandomScenarios < matRad_RandomScenariosAbstract
    %  matRad_RandomScenarios
    %  Implements randomly sampled scenarios.
    %
    % constructor
    %   matRad_RandomScenarios()
    %   matRad_RandomScenarios(ct)
    %
    % input:
    %   ct:                 ct cube
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2022-2026 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (SetAccess = protected)
        name = 'Randomly sampled Scenarios'
        shortName   = 'rndScen'
    end

    methods

        function this = matRad_RandomScenarios(ct)
            if nargin == 0
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end

            this@matRad_RandomScenariosAbstract(superclassArgs{:});
        end

    end

end
