classdef matRad_TruncatedImportanceScenarios < matRad_GriddedScenariosAbstract
    %  matRad_TruncatedImportanceScenarios
    %  Implements gridded importance scenarios truncated to the wcSigma
    %  normalized uncertainty radius. The grid is first created like
    %  matRad_ImportanceScenarios and then scenarios outside the configured
    %  normalized radius are removed.
    %
    % constructor
    %   matRad_TruncatedImportanceScenarios()
    %   matRad_TruncatedImportanceScenarios(ct)
    %
    % input:
    %   ct:                 ct cube
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
        numOfSetupGridPoints = 9
        numOfRangeGridPoints = 9
    end

    properties (SetAccess = protected)
        shortName   = 'truncatedImpScen'
        name        = 'Truncated Gridded Scenarios with Importance Weights'
    end

    methods

        function this = matRad_TruncatedImportanceScenarios(ct)
            if nargin == 0
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end

            this@matRad_GriddedScenariosAbstract(superclassArgs{:});

            this.updateScenarios();
        end

        function set.numOfSetupGridPoints(this, numGridPoints)
            valid = isscalar(numGridPoints) && numGridPoints > 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid number of setup grid points, needs to be a positive scalar!');
            end
            this.numOfSetupGridPoints = numGridPoints;
            this.updateScenarios();
        end

        function set.numOfRangeGridPoints(this, numGridPoints)
            valid = isscalar(numGridPoints) && numGridPoints > 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid number of range grid points, needs to be a positive scalar!');
            end
            this.numOfRangeGridPoints = numGridPoints;
            this.updateScenarios();
        end

        function scenarios = updateScenarios(this)
            updateScenarios@matRad_GriddedScenariosAbstract(this);

            normalizedRadius = this.normalizedScenarioRadius();
            tolerance = 100 * eps(max(1, this.wcSigma));
            keepIx = normalizedRadius <= this.wcSigma + tolerance;
            keepIx = keepIx | all(abs(this.scenarioValues) <= eps, 2);

            if ~all(keepIx)
                this.applyScenarioFilter(keepIx);
            end

            scenarios = this.scenForProb;
        end

    end

    methods (Access = private)

        function normalizedRadius = normalizedScenarioRadius(this)
            scenarioValues = this.scenarioValues;
            scenarioScale = [this.scenarioComponents.scale];
            activeIx = [this.scenarioComponents.active];
            normalizedRadius = zeros(size(scenarioValues, 1), 1);

            if any(activeIx)
                activeValues = scenarioValues(:, activeIx);
                activeScale = scenarioScale(activeIx);
                scaledValues = bsxfun(@rdivide, activeValues, activeScale);
                normalizedRadius = sqrt(sum(scaledValues.^2, 2));
            end
        end

        function applyScenarioFilter(this, keepIx)
            scenForProb = this.scenForProb(keepIx, :);
            linearMask = this.linearMask(keepIx, :);
            scenarioValues = this.scenarioValues(keepIx, :);
            ctScenIds = scenForProb(:, 1);
            scenProb = this.scenProb(keepIx);
            scenWeight = scenProb ./ sum(scenProb);

            scenMaskSize = [this.numOfAvailableCtScen, ...
                            this.totNumShiftScen, ...
                            this.totNumRangeScen];
            scenMask = false(scenMaskSize);
            maskIx = sub2ind(scenMaskSize, linearMask(:, 1), ...
                             linearMask(:, 2), linearMask(:, 3));
            scenMask(maskIx) = true;

            this.setScenarioRealizations(this.scenarioComponents, scenarioValues, ...
                                         ctScenIds, scenProb, scenWeight, scenForProb, ...
                                         linearMask, scenMaskSize, 'legacy-grid');
        end

    end

end
