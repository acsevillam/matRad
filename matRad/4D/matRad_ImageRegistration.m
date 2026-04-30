classdef (Abstract) matRad_ImageRegistration
    % matRad_ImageRegistration Abstract superclass for image registration
    %   Defines the common interface for image registration classes that
    %   compute deformation vector fields and propagate contours between CT
    %   scenarios.
    %
    % References
    %   -
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


    properties (Abstract, Constant)
        name                %Display name of the registration method. Needs to be implemented in sub-classes.
    end

    properties (Abstract)
        ct
        cst
        refScen
        metadata
    end

    methods
        function obj = matRad_ImageRegistration(dataStruct)
            if nargin > 0 && ~isempty(dataStruct) && isstruct(dataStruct)
                obj = assignCommonPropertiesFromStruct(obj,dataStruct);
            end
        end

        % Overload the struct function to return the common registration data.
        function s = struct(obj)
            s.className = class(obj);
            s.ct = obj.ct;
            s.cst = obj.cst;
            s.refScen = obj.refScen;
            s.metadata = obj.metadata;
        end
    end

    methods (Access = private)
        function obj = assignCommonPropertiesFromStruct(obj,s)
            for fn = fieldnames(s)'
                try
                    obj.(fn{1}) = s.(fn{1});
                catch
                    continue;
                end
            end
        end
    end

    methods (Abstract)
        ct = calcDVF(obj)
        [ct,cst] = propContours(obj)
    end

end
