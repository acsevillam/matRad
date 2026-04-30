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

            %TODO: We could do this automatically in the superclass
            %Octave 5 has a bug there and throws an error
            this.updateScenarios();
        end

        function set.numOfSetupGridPoints(this,numGridPoints)
            valid = isscalar(numGridPoints) && numGridPoints > 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid number of setup grid points, needs to be a positive scalar!');
            end
            this.numOfSetupGridPoints = numGridPoints;
            this.updateScenarios();
        end

        function set.numOfRangeGridPoints(this,numGridPoints)
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

            truncationRadius = this.wcSigma;

            scenValues = this.scenForProb(:,2:6);
            scenScale = [this.shiftSD this.rangeAbsSD this.rangeRelSD./100];
            scenScale(scenScale == 0) = eps;

            normalizedRadius = sqrt(sum(bsxfun(@rdivide,scenValues,scenScale).^2,2));
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
            this.scenForProb = this.scenForProb(keepIx,:);
            this.linearMask = this.linearMask(keepIx,:);
            this.ctScenIx = this.scenForProb(:,1);
            this.totNumScen = size(this.scenForProb,1);

            this.relRangeShift = this.scenForProb(:,6);
            this.absRangeShift = this.scenForProb(:,5);
            this.isoShift = this.scenForProb(:,2:4);
            this.maxAbsRangeShift = max(abs(this.absRangeShift));
            this.maxRelRangeShift = max(abs(this.relRangeShift));

            this.scenProb = this.scenProb(keepIx);
            this.scenWeight = this.scenProb./sum(this.scenProb);

            this.scenMask = false(this.numOfAvailableCtScen, ...
                this.totNumShiftScen,this.totNumRangeScen);
            maskIx = sub2ind(size(this.scenMask), ...
                this.linearMask(:,1),this.linearMask(:,2),this.linearMask(:,3));
            this.scenMask(maskIx) = true;
        end
    end
end
