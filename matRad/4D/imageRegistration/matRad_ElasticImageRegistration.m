classdef matRad_ElasticImageRegistration < matRad_ImageRegistration
    % matRad_ElasticImageRegistration Demons-based elastic image registration
    %
    % call
    %   obj = matRad_ElasticImageRegistration(ct,cst)
    %   obj = matRad_ElasticImageRegistration(ct,cst,refScen)
    %   obj = matRad_ElasticImageRegistration(ct,cst,refScen,metadata)
%   obj = matRad_ElasticImageRegistration(dataStruct)
%
%   This registration method requires MATLAB with the Image Processing
%   Toolbox. It is not available in Octave.
    %
    % input
    %   ct:             matRad ct struct
    %   cst:            matRad cst struct
    %   refScen:        reference CT scenario index (default: 1)
    %   metadata:       struct with optional registration settings
    %                   metadata.dvfType:      'pull' or 'push' (default: 'pull')
    %                   metadata.nItera:       number of registration iterations (default: 100)
    %                   metadata.pyramLevels:  number of pyramid levels (default: 1)
    %                   metadata.smoothLevels: accumulated field smoothing (default: 1)
    %                   metadata.dvfUnits:     DVF units, currently 'voxel'
    %
    % output
    %   obj:            elastic image registration object
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


    properties (Constant)
        name = 'Elastic Registration';
    end

    properties
        ct
        cst
        refScen
        metadata
    end

    methods (Access = public)


        function obj = matRad_ElasticImageRegistration(ct,cst,refScen,metadata)

            if nargin == 1 && isstruct(ct)
                inputStruct = ct;
                initFromStruct = true;
            else
                inputStruct = [];
                initFromStruct = false;
            end

            obj@matRad_ImageRegistration(inputStruct);

            if initFromStruct
                return;
            end

            if nargin < 2
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('ct and cst must be provided for elastic image registration.');
            end

            if nargin < 3 || isempty(refScen)
                refScen = 1;
            end

            refScen = matRad_ElasticImageRegistration.validateReferenceScenario(ct,refScen);

            if nargin < 4 || isempty(metadata)
                metadata = struct();
            end

            metadata = matRad_ElasticImageRegistration.setDefaultMetadata(metadata);

            obj.ct = ct;
            obj.cst = cst;
            obj.refScen = refScen;
            obj.metadata = metadata;

        end

        % Calculate the deformation vector fields.
        %
        % output
        %   ct: matRad ct struct with deformation vector fields
        function ct = calcDVF(obj)

            % Non rigid registration demons-based. Calculates the DVF(Displacement
            % Vector Field) that models the transformation.
            matRad_cfg = MatRad_Config.instance();
            matRad_ElasticImageRegistration.checkImageProcessingDependency(false);

            nScen = obj.ct.numOfCtScen;
            obj.ct.dvf = cell(1,nScen);
            obj.ct.dvfType = obj.metadata.dvfType;
            obj.ct.dvfUnits = obj.metadata.dvfUnits;
            obj.ct.dvfDim = obj.ct.cubeDim;
            obj.ct.refScen = obj.refScen;
            obj.ct.dvfMetadata = obj.metadata;
            obj.ct.dvfMetadata.refScen = obj.refScen;
            obj.ct.dvfMetadata.referenceCtScen = obj.refScen;

            switch obj.metadata.dvfType
                case 'push'
                    for scen = 1:nScen
                        matRad_cfg.dispInfo('Registering scenario %d.\n',scen);
                        [obj.ct.dvf{scen},~] = imregdemons(obj.ct.cubeHU{obj.refScen}, ...
                            obj.ct.cubeHU{scen},obj.metadata.nItera, ...
                            'PyramidLevels',obj.metadata.pyramLevels, ...
                            'AccumulatedFieldSmoothing',obj.metadata.smoothLevels);
                        obj.ct.dvf{scen} = permute(obj.ct.dvf{scen},[4 1 2 3]);
                    end
                case 'pull'
                    for scen = 1:nScen
                        matRad_cfg.dispInfo('Registering scenario %d.\n',scen);
                        [obj.ct.dvf{scen},~] = imregdemons(obj.ct.cubeHU{scen}, ...
                            obj.ct.cubeHU{obj.refScen},obj.metadata.nItera, ...
                            'PyramidLevels',obj.metadata.pyramLevels, ...
                            'AccumulatedFieldSmoothing',obj.metadata.smoothLevels);
                        obj.ct.dvf{scen} = permute(obj.ct.dvf{scen},[4 1 2 3]);
                    end
            end

            ct = obj.ct;

        end

        % Calculate the deformation vector fields on a resized CT grid.
        %
        % input
        %   resolution: struct with x, y, and z grid spacing
        %
        % output
        %   ct:         matRad ct struct with deformation vector fields
        function ct = calcDVF_resized(obj,resolution)

            matRad_cfg = MatRad_Config.instance();
            matRad_ElasticImageRegistration.checkImageProcessingDependency(false);
            nScen = obj.ct.numOfCtScen;
            cubeHU_resized = cell(1,nScen);
            tmpGrid.x = obj.ct.x(1):resolution.x:obj.ct.x(end);
            tmpGrid.y = obj.ct.y(1):resolution.y:obj.ct.y(end);
            tmpGrid.z = obj.ct.z(1):resolution.z:obj.ct.z(end);
            obj.ct.dvf = cell(1,nScen);
            obj.ct.dvfType = obj.metadata.dvfType;
            obj.ct.dvfUnits = obj.metadata.dvfUnits;
            obj.ct.dvfDim = [numel(tmpGrid.y) numel(tmpGrid.x) numel(tmpGrid.z)];
            obj.ct.refScen = obj.refScen;
            obj.ct.dvfMetadata = obj.metadata;
            obj.ct.dvfMetadata.refScen = obj.refScen;
            obj.ct.dvfMetadata.referenceCtScen = obj.refScen;

            for scen = 1:nScen
                cubeHU_resized{scen} = matRad_interp3(obj.ct.x,obj.ct.y,obj.ct.z, ...
                                       obj.ct.cubeHU{scen}, ...
                                       tmpGrid.x,tmpGrid.y',tmpGrid.z,'nearest');
            end

            % Non rigid registration demons-based. Calculates the DVF(Displacement
            % Vector Field) that models the transformation.

            switch obj.metadata.dvfType
                case 'push'
                    for scen = 1:nScen
                        matRad_cfg.dispInfo('Registering scenario %d.\n',scen);
                        [obj.ct.dvf{scen},~] = imregdemons(cubeHU_resized{obj.refScen}, ...
                            cubeHU_resized{scen},obj.metadata.nItera, ...
                            'PyramidLevels',obj.metadata.pyramLevels, ...
                            'AccumulatedFieldSmoothing',obj.metadata.smoothLevels);
                        obj.ct.dvf{scen} = permute(obj.ct.dvf{scen},[4 1 2 3]);
                    end
                case 'pull'
                    for scen = 1:nScen
                        matRad_cfg.dispInfo('Registering scenario %d.\n',scen);
                        [obj.ct.dvf{scen},~] = imregdemons(cubeHU_resized{scen}, ...
                            cubeHU_resized{obj.refScen},obj.metadata.nItera, ...
                            'PyramidLevels',obj.metadata.pyramLevels, ...
                            'AccumulatedFieldSmoothing',obj.metadata.smoothLevels);
                        obj.ct.dvf{scen} = permute(obj.ct.dvf{scen},[4 1 2 3]);
                    end
            end

            ct = obj.ct;

        end

        % Propagate contours using push deformation vector fields.
        %
        % output
        %   ct:  matRad ct struct with deformation vector fields
        %   cst: matRad cst struct with propagated contours
        function [ct,cst] = propContours(obj)

            matRad_cfg = MatRad_Config.instance();
            matRad_ElasticImageRegistration.checkImageProcessingDependency(true);

            if ~strcmp(obj.metadata.dvfType,'push')
                matRad_cfg.dispError('Contour propagation requires push DVFs.');
            end

            if ~isfield(obj.ct,'dvf') || ~isfield(obj.ct,'dvfType')
                obj.ct = obj.calcDVF();
            end

            [numOfStruct, ~] = size(obj.cst);
            for structure = 1:numOfStruct

                if numel(obj.cst{structure,4}) < obj.refScen
                    matRad_cfg.dispError('Structure %d does not contain contours for reference CT scenario %d.', ...
                        structure,obj.refScen);
                end

                if ~isempty(obj.cst{structure,4}{1,obj.refScen})

                    % Obtaining the fixed cubic structure from the linear indices
                    cubeHU_fixed = zeros(obj.ct.cubeDim);
                    cst_fixed = obj.cst{structure,4}{1,obj.refScen};
                    [x,y,z] = ind2sub(obj.ct.cubeDim,cst_fixed);

                    % The HU value of the corresponding tomography is assigned to each position
                    for j = 1:length(x)
                        cubeHU_fixed(x(j),y(j),z(j)) = 1; %obj.cubeHU{1}(x(j),y(j),z(j));
                    end

                    % The DVF transformation is applied and the linear values are found
                    matRad_cfg.dispInfo('Propagating contours of structure %d.\n',structure);
                    for scen = 1:obj.ct.numOfCtScen
                        if scen == obj.refScen
                            continue;
                        end
                        cubeHU_estimated = imwarp(cubeHU_fixed, permute(obj.ct.dvf{scen},[2 3 4 1]));
                        obj.cst{structure,4}{1,scen} = find(cubeHU_estimated);
                    end
                else
                    matRad_cfg.dispWarning('Structure %d has no reference contour. Leaving propagated contours empty.\n',structure);
                    for scen = 1:obj.ct.numOfCtScen
                        if scen ~= obj.refScen
                            obj.cst{structure,4}{1,scen} = [];
                        end
                    end
                end
            end

            ct = obj.ct;
            cst = obj.cst;

        end

    end

    methods (Static, Access = private)
        function checkImageProcessingDependency(requireWarp)

            matRad_cfg = MatRad_Config.instance();
            env = matRad_getEnvironment();

            if strcmp(env,'OCTAVE')
                matRad_cfg.dispError(['matRad_ElasticImageRegistration requires MATLAB ', ...
                    'with the Image Processing Toolbox.']);
            end

            if exist('imregdemons','file') ~= 2
                matRad_cfg.dispError(['matRad_ElasticImageRegistration requires ', ...
                    'imregdemons from the Image Processing Toolbox.']);
            end

            if requireWarp && exist('imwarp','file') ~= 2
                matRad_cfg.dispError(['Contour propagation requires imwarp from ', ...
                    'the Image Processing Toolbox.']);
            end

            if exist('license','builtin') == 5 || exist('license','file') == 2
                hasImageToolboxLicense = true;
                try
                    hasImageToolboxLicense = license('test','image_toolbox');
                catch
                end
                if ~hasImageToolboxLicense
                    matRad_cfg.dispError(['matRad_ElasticImageRegistration requires ', ...
                        'an available Image Processing Toolbox license.']);
                end
            end
        end

        function refScen = validateReferenceScenario(ct,refScen)

            matRad_cfg = MatRad_Config.instance();

            if ~(isnumeric(refScen) && isscalar(refScen) && isfinite(refScen) && ...
                    round(refScen) == refScen && refScen >= 1)
                matRad_cfg.dispError('refScen must be a positive integer scalar.');
            end

            if isfield(ct,'numOfCtScen') && refScen > ct.numOfCtScen
                matRad_cfg.dispError('refScen (%d) exceeds ct.numOfCtScen (%d).', ...
                    refScen,ct.numOfCtScen);
            end

            refScen = double(refScen);
        end

        function metadata = setDefaultMetadata(metadata)

            matRad_cfg = MatRad_Config.instance();

            if ~isstruct(metadata)
                matRad_cfg.dispError('metadata must be a struct.');
            end

            if ~isfield(metadata,'dvfType') || isempty(metadata.dvfType)
                metadata.dvfType = 'pull';
            else
                metadata.dvfType = char(metadata.dvfType);
            end

            if ~any(strcmp(metadata.dvfType,{'pull','push'}))
                matRad_cfg.dispError('metadata.dvfType must be ''pull'' or ''push''.');
            end

            if ~isfield(metadata,'nItera') || isempty(metadata.nItera)
                metadata.nItera = 100;
            end

            if ~isfield(metadata,'pyramLevels') || isempty(metadata.pyramLevels)
                metadata.pyramLevels = 1;
            end

            if ~isfield(metadata,'smoothLevels') || isempty(metadata.smoothLevels)
                metadata.smoothLevels = 1;
            end

            if ~isfield(metadata,'dvfUnits') || isempty(metadata.dvfUnits)
                metadata.dvfUnits = 'voxel';
            else
                metadata.dvfUnits = char(metadata.dvfUnits);
            end

            if ~any(strcmp(metadata.dvfUnits,{'voxel','mm'}))
                matRad_cfg.dispError('metadata.dvfUnits must be ''voxel'' or ''mm''.');
            end
        end
    end

end
