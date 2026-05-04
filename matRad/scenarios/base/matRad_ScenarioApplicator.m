classdef matRad_ScenarioApplicator < handle
% matRad_ScenarioApplicator base class for scenario realization applicators.
%
% Applicators separate the scenario values from the physical place where a
% realization modifies dose calculation.
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

    properties (SetAccess = protected)
        name = 'base';
        componentNames = {};
    end

    methods
        function tf = supportsComponent(this,componentName)
            tf = any(strcmp(this.componentNames,componentName));
        end
    end
end
