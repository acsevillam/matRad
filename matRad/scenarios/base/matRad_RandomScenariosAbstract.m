classdef (Abstract) matRad_RandomScenariosAbstract < matRad_ScenarioModel
    %  matRad_RandomScenariosAbstract
    %  Shared implementation for randomly sampled scenario models.
    %
    % constructor
    %   matRad_RandomScenariosAbstract()
    %   matRad_RandomScenariosAbstract(ct)
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

    properties
        includeNominalScenario = false  % Forces inclusion of the nominal scenario
        nSamples = 10                   % Standard number of random samples
        randomSeed = []                 % optional seed for reproducible random scenario values
    end

    % Deprecated Properties that were used
    properties (Dependent)
        numOfShiftScen
        numOfRangeShiftScen
    end

    methods

        function this = matRad_RandomScenariosAbstract(ct)
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

        %% Setters & Update
        function set.nSamples(this, nSamples)
            valid = isnumeric(nSamples) && isscalar(nSamples) && mod(nSamples, 1) == 0 && nSamples > 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for nSamples! Needs to be a positive integer!');
            end
            this.nSamples = nSamples;
            this.updateScenarios();
        end

        function set.randomSeed(this, randomSeed)
            valid = isempty(randomSeed) || (isnumeric(randomSeed) && isscalar(randomSeed) && ...
                                            isfinite(randomSeed) && round(randomSeed) == randomSeed && randomSeed >= 0);
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for randomSeed! Needs to be empty or a non-negative integer scalar!');
            end
            this.randomSeed = randomSeed;
            this.updateScenarios();
        end

        function set.numOfShiftScen(this, numOfShiftScen)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfShiftScen of the scenario class will soon be deprecated! Use nSamples instead');

            % That's for downwards compatibility
            if ~isscalar(numOfShiftScen)
                numOfShiftScen = unique(numOfShiftScen);
            end

            this.nSamples = numOfShiftScen;
        end

        function value = get.numOfShiftScen(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfShiftScen of the scenario class will soon be deprecated! Use nSamples instead');
            value = this.nSamples;
        end

        function set.numOfRangeShiftScen(this, numOfRangeShiftScen)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfRangeShiftScen of the scenario class will soon be deprecated! Use nSamples instead');
            this.nSamples = numOfRangeShiftScen;
        end

        function value = get.numOfRangeShiftScen(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfRangeShiftScen of the scenario class will soon be deprecated! Use nSamples instead');
            value = this.nSamples;
        end

        function scenarios = updateScenarios(this)
            matRad_cfg = MatRad_Config.instance();

            this.numOfCtScen = size(this.ctScenProb, 1);
            components = this.getScenarioComponents();
            activeIx = [components.active];
            nScenarioSamples = helper_scenarioCountForActiveDimension(any(activeIx), this.nSamples);

            scenarios = this.sampleScenarioValues(nScenarioSamples, ...
                                                  components, activeIx);

            if this.includeNominalScenario
                % We include the nominal scenario by just replacing the
                % first one to keep the number of scenarios the same.
                scenarios(1, :) = 0;
            end

            % Handle 4D CT scenario ids.
            ctScenIds = repmat(this.ctScenProb(:, 1)', size(scenarios, 1), 1);
            ctScenIds = ctScenIds(:);
            scenarioValues = repmat(scenarios, [this.numOfCtScen 1]);
            scenarios = horzcat(ctScenIds, scenarioValues);

            scenProb = matRad_computeScenarioProbabilities(components, scenarioValues, ...
                                                           this.ctScenProb, ctScenIds);

            % Scenario weight.
            scenWeight = [];
            for sCt = 1:this.numOfCtScen
                % Equal weights within a phase since they have been randomly
                % sampled; this is approximate when the nominal scenario is forced.
                scenWeight = [scenWeight; ...
                              this.ctScenProb(sCt, 2) * ones(nScenarioSamples, 1) ./ nScenarioSamples];
            end

            setupActive = any(strcmp({components.applicator}, 'setup') & [components.active]);
            rangeActive = any(strcmp({components.applicator}, 'range') & [components.active]);
            gantryActive = any(strcmp({components.applicator}, 'gantry') & [components.active]);
            couchActive = any(strcmp({components.applicator}, 'couch') & [components.active]);
            this.totNumShiftScen = helper_scenarioCountForActiveDimension(setupActive, this.nSamples);
            this.totNumRangeScen = helper_scenarioCountForActiveDimension(rangeActive, this.nSamples);
            this.totNumGantryScen = helper_scenarioCountForActiveDimension(gantryActive, this.nSamples);
            this.totNumCouchScen = helper_scenarioCountForActiveDimension(couchActive, this.nSamples);
            this.totNumScen = nScenarioSamples * this.numOfCtScen;

            linearMaskAll = zeros(this.totNumScen, 5);
            rowIx = 0;
            for sCt = 1:this.numOfCtScen
                ctScenId = this.ctScenProb(sCt, 1);
                for sampleIx = 1:nScenarioSamples
                    rowIx = rowIx + 1;
                    linearMaskAll(rowIx, :) = [ctScenId, ...
                                               helper_scenarioIndexForSample(setupActive, sampleIx), ...
                                               helper_scenarioIndexForSample(rangeActive, sampleIx), ...
                                               helper_scenarioIndexForSample(gantryActive, sampleIx), ...
                                               helper_scenarioIndexForSample(couchActive, sampleIx)];
                end
            end

            if this.hasOnlyLegacyScenarioDimensions()
                storagePolicy = 'legacy-grid';
                storageSize = [this.numOfAvailableCtScen, this.totNumShiftScen, ...
                               this.totNumRangeScen];
                storageSubscripts = linearMaskAll(:, 1:3);
            else
                storagePolicy = 'compact-realization';
                storageSize = [this.totNumScen 1];
                storageSubscripts = [(1:this.totNumScen)' ones(this.totNumScen, 1)];
            end

            this.setScenarioRealizations(components, scenarioValues, ctScenIds, ...
                                         scenProb, scenWeight, scenarios, storageSubscripts, storageSize, ...
                                         storagePolicy);
            totNumScen = sum(this.scenMask(:));

            if totNumScen ~= this.totNumScen
                matRad_cfg.dispWarning('Check Implementation of Total Scenario computation - given %d but found %d!', this.totNumScen, totNumScen);
                this.totNumScen = totNumScen;
            end
        end

    end

    methods (Access = protected)

        function scenarios = sampleScenarioValues(this, nScenarioSamples, components, activeIx)
            % Multivariate normal sampling over active components only.
            scenarios = zeros(nScenarioSamples, numel(components));
            if any(activeIx)
                scales = [components(activeIx).scale];
                randomValues = this.sampleStandardNormalValues(nScenarioSamples, ...
                                                               sum(activeIx));
                scenarios(:, activeIx) = bsxfun(@times, randomValues, scales);
            end
        end

        function randomValues = sampleStandardNormalValues(this, numRows, numColumns)
            if isempty(this.randomSeed)
                randomValues = randn(numRows, numColumns);
            else
                rngState = rng;
                cleanupObj = onCleanup(@() rng(rngState));
                rng(this.randomSeed);
                randomValues = randn(numRows, numColumns);
                clear cleanupObj;
            end
        end

        function validateScenarioDimensionSupport(~, ~)
            % Random sampling supports all currently declared continuous
            % scenario dimensions.
        end

    end

end

function n = helper_scenarioCountForActiveDimension(isActive, nSamples)

    if isActive
        n = nSamples;
    else
        n = 1;
    end

end

function scenIx = helper_scenarioIndexForSample(isActive, sampleIx)

    if isActive
        scenIx = sampleIx;
    else
        scenIx = 1;
    end

end
