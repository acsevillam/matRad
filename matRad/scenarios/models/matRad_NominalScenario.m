classdef matRad_NominalScenario < matRad_ScenarioModel
    %  matRad_RandomScenarios
    %  Implements a single nominal planning scenario
    %
    % constructor
    %   matRad_NominalScenario()
    %   matRad_NominalScenario(ct)
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
        shortName   = 'nomScen'
        name        = 'Nominal Scenario'
    end

    methods

        function this = matRad_NominalScenario(ct)
            if nargin == 0
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end
            this@matRad_ScenarioModel(superclassArgs{:});

            % TODO: We could do this automatically in the superclass
            % Octave 5 has a bug there and throws an error
            this.updateScenarios();
        end

        function scenarios = updateScenarios(this)
            this.numOfCtScen = size(this.ctScenProb, 1);

            % set variables
            this.totNumShiftScen = 1;
            this.totNumRangeScen = 1;
            this.totNumGantryScen = 1;
            this.totNumCouchScen = 1;
            this.totNumScen = this.numOfCtScen;

            components = this.getScenarioComponents();
            ctScenIds = this.ctScenProb(:, 1);
            scenarioValues = zeros(this.numOfCtScen, numel(components));
            scenarios = [ctScenIds scenarioValues];

            % Get Scenario probability
            scenProb = matRad_computeScenarioProbabilities(components, scenarioValues, ...
                                                           this.ctScenProb, ctScenIds);

            % Get relative (normalized) weight of the scenario
            scenWeight = scenProb ./ sum(scenProb);

            if this.hasOnlyLegacyScenarioDimensions()
                storagePolicy = 'legacy-grid';
                storageSize = [this.numOfAvailableCtScen, this.totNumShiftScen, this.totNumRangeScen];
                storageSubscripts = [ctScenIds ones(this.numOfCtScen, 2)];
            else
                storagePolicy = 'compact-realization';
                storageSize = [this.numOfCtScen 1];
                storageSubscripts = [(1:this.numOfCtScen)' ones(this.numOfCtScen, 1)];
            end

            this.setScenarioRealizations(components, scenarioValues, ctScenIds, ...
                                         scenProb, scenWeight, scenarios, storageSubscripts, storageSize, ...
                                         storagePolicy);

            if sum(this.scenMask(:)) ~= this.totNumScen
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning(['Check Implementation of Total Scenario computation - given %d ' ...
                                        'but found %d!'], this.totNumScen, sum(this.scenMask(:)));
                this.totNumScen = sum(this.scenMask(:));
            end
        end

    end

    methods (Access = protected)

        function validateScenarioDimensionSupport(~, ~)
            % Nominal scenarios can represent any registered dimension at its
            % nominal value.
        end

    end
end
