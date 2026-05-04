classdef matRad_TruncatedImportanceScenarios < matRad_GriddedScenariosAbstract
%  matRad_TruncatedImportanceScenarios
%  Implements gridded importance scenarios truncated to the wcSigma
%  Mahalanobis ball. The grid is first created like matRad_ImportanceScenarios
%  and then scenarios outside the configured uncertainty radius are removed.
%  This is useful for combined setup/range grids where corner scenarios can
%  otherwise exceed the intended robust uncertainty radius.
%
% constructor
%   matRad_TruncatedImportanceScenarios()
%   matRad_TruncatedImportanceScenarios(ct)
%
% input
%   ct:                 ct cube
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

    properties (AbortSet = true)
        numOfSetupGridPoints = 9;
        numOfRangeGridPoints = 9;
    end

    properties (SetAccess=protected)
        name = 'truncatedImpScen';
    end

    methods
        function this = matRad_TruncatedImportanceScenarios(ct)
            if nargin == 0
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end

            this@matRad_GriddedScenariosAbstract(superclassArgs{:});

            this.initializeScenarioModel();
        end

        function set.numOfSetupGridPoints(this,numGridPoints)
            valid = isscalar(numGridPoints) && numGridPoints > 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid number of setup grid points, needs to be a positive scalar!');
            end
            this.numOfSetupGridPoints = numGridPoints;
            this.requestScenarioUpdate();
        end

        function set.numOfRangeGridPoints(this,numGridPoints)
            valid = isscalar(numGridPoints) && numGridPoints > 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid number of range grid points, needs to be a positive scalar!');
            end
            this.numOfRangeGridPoints = numGridPoints;
            this.requestScenarioUpdate();
        end

        function scenarios = updateScenarios(this)
            updateScenarios@matRad_GriddedScenariosAbstract(this);

            truncationRadius = this.wcSigma;

            scenValues = this.scenarioValues;
            scenScale = [this.scenarioComponents.scale];
            activeIx = [this.scenarioComponents.active];
            normalizedRadius = zeros(size(scenValues,1),1);
            if any(activeIx)
                normalizedRadius = sqrt(sum(bsxfun(@rdivide, ...
                    scenValues(:,activeIx),scenScale(activeIx)).^2,2));
            end
            keepIx = normalizedRadius <= truncationRadius + 100*eps(max(1,truncationRadius));
            keepIx = keepIx | all(abs(scenValues) <= eps,2);

            if ~all(keepIx)
                this.applyScenarioFilter(keepIx);
            end

            scenarios = this.scenForProb;
        end
    end

    methods (Access = private)
        function applyScenarioFilter(this,keepIx)
            scenForProb = this.scenForProb(keepIx,:);
            linearMask = this.linearMask(keepIx,:);
            scenarioValues = this.scenarioValues(keepIx,:);
            ctScenIds = scenForProb(:,1);
            this.totNumScen = size(scenForProb,1);

            scenProb = this.scenProb(keepIx);
            scenWeight = this.scenWeight(keepIx);

            scenMask = false(this.numOfAvailableCtScen, ...
                this.totNumShiftScen,this.totNumRangeScen);
            maskIx = sub2ind(size(scenMask), ...
                linearMask(:,1),linearMask(:,2),linearMask(:,3));
            scenMask(maskIx) = true;
            [scenarioValues,ctScenIds,scenProb,scenWeight,scenForProb,linearMask,scenMask] = ...
                matRad_filterZeroProbabilityScenarios(scenarioValues,ctScenIds, ...
                scenProb,scenWeight,scenForProb,linearMask,scenMask);
            this.totNumScen = size(scenarioValues,1);
            this.setScenarioRealizations(this.scenarioComponents,scenarioValues,ctScenIds, ...
                scenProb,scenWeight,scenForProb,linearMask,scenMask);
        end
    end
end
