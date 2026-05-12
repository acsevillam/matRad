classdef matRad_RandomScenarios < matRad_ScenarioModel
%  matRad_RandomScenarios
%  Implements randomly sampled scenarios
%
% constructor
%   matRad_RandomScenarios()
%   matRad_RandomScenarios(ct)
%
% input
%   ct:                 ct cube
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
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
        includeNominalScenario = false; %Forces inclusion of the nominal scenario
        nSamples = 10;                  %Standard number of random samples
        randomSeed = [];                %Optional seed for reproducible random scenario values
    end

    %Deprecated Properties that were used
    properties (Dependent)
        numOfShiftScen;
        numOfRangeShiftScen;
    end

    properties (SetAccess=protected)
        name = 'rndScen';
    end

    methods
        function this = matRad_RandomScenarios(ct)
            if nargin == 0
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end

            this@matRad_ScenarioModel(superclassArgs{:});

            this.initializeScenarioModel();
        end

        %% Setters & Update
        function set.includeNominalScenario(this,includeNominalScenario)
            valid = isscalar(includeNominalScenario) && ...
                (islogical(includeNominalScenario) || isnumeric(includeNominalScenario));
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for includeNominalScenario! Needs to be a boolean / logical value!');
            end
            this.includeNominalScenario = includeNominalScenario;
            this.requestScenarioUpdate();
        end

        function set.nSamples(this,nSamples)
            valid = isnumeric(nSamples) && isscalar(nSamples) && mod(nSamples,1) == 0 && nSamples > 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for nSamples! Needs to be a positive integer!');
            end
            this.nSamples = nSamples;
            this.requestScenarioUpdate();
        end

        function set.randomSeed(this,randomSeed)
            valid = isempty(randomSeed) || (isnumeric(randomSeed) && isscalar(randomSeed) && ...
                isfinite(randomSeed) && round(randomSeed) == randomSeed && randomSeed >= 0);
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for randomSeed! Needs to be empty or a non-negative integer scalar!');
            end
            this.randomSeed = randomSeed;
            this.requestScenarioUpdate();
        end

        function set.numOfShiftScen(this,numOfShiftScen)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfShiftScen of the scenario class will soon be deprecated! Use nSamples instead');

            %That's for downwards compatibility
            if ~isscalar(numOfShiftScen)
                numOfShiftScen = unique(numOfShiftScen);
            end

            this.nSamples = numOfShiftScen;
        end

        function  value = get.numOfShiftScen(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfShiftScen of the scenario class will soon be deprecated! Use nSamples instead');
            value = this.nSamples;
        end

        function set.numOfRangeShiftScen(this,numOfRangeShiftScen)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfRangeShiftScen of the scenario class will soon be deprecated! Use nSamples instead');
            this.nSamples = numOfRangeShiftScen;
        end

        function  value = get.numOfRangeShiftScen(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfRangeShiftScen of the scenario class will soon be deprecated! Use nSamples instead');
            value = this.nSamples;
        end

        function scenarios = updateScenarios(this)
            matRad_cfg = MatRad_Config.instance();

            this.numOfCtScen = size(this.ctScenProb,1);
            components = this.getScenarioComponents();
            activeIx = [components.active];
            nScenarioSamples = scenarioCountForActiveDimension(any(activeIx),this.nSamples);

            %Multivariate Normal Sampling over active components only.
            scenarios = zeros(nScenarioSamples,numel(components));
            if any(activeIx)
                scales = [components(activeIx).scale];
                Sigma = diag(scales.^2);
                cs = chol(Sigma);
                if isempty(this.randomSeed)
                    randomValues = randn(nScenarioSamples,sum(activeIx));
                else
                    rngState = rng;
                    cleanupObj = onCleanup(@() rng(rngState));
                    rng(this.randomSeed);
                    randomValues = randn(nScenarioSamples,sum(activeIx));
                    clear cleanupObj;
                end
                scenarios(:,activeIx) = randomValues * cs;
            end

            if this.includeNominalScenario
                %We include the nominal scenario by just replacing the
                %first one to keep the number of scenarios the same
                scenarios(1,:) = 0;
            end

            %Handle 4D CT scenario ids
            ctScenIds = repmat(this.ctScenProb(:,1)',size(scenarios,1),1);
            ctScenIds = ctScenIds(:);
            scenarios = horzcat(ctScenIds, repmat(scenarios,[this.numOfCtScen 1]));

            %Scenario weight
            scenWeight = [];
            for sCt = 1:this.numOfCtScen
                %equal weights within a phase since they have been randomly sampled
                %(not entirely true if the Nominal scenario was forced)
                scenWeight = [scenWeight; this.ctScenProb(sCt,2) * ones(nScenarioSamples,1)./nScenarioSamples];
            end

            %set variables
            setupActive = any(strcmp({components.applicator},'setup') & [components.active]);
            rangeActive = any(strcmp({components.applicator},'range') & [components.active]);
            gantryActive = any(strcmp({components.applicator},'gantry') & [components.active]);
            couchActive = any(strcmp({components.applicator},'couch') & [components.active]);
            this.totNumShiftScen = scenarioCountForActiveDimension(setupActive,this.nSamples);
            this.totNumRangeScen = scenarioCountForActiveDimension(rangeActive,this.nSamples);
            this.totNumGantryScen = scenarioCountForActiveDimension(gantryActive,this.nSamples);
            this.totNumCouchScen = scenarioCountForActiveDimension(couchActive,this.nSamples);
            this.totNumScen = nScenarioSamples * this.numOfCtScen; %check because of CT scenarios

            linearMaskAll = zeros(this.totNumScen,5);
            rowIx = 0;
            for sCt = 1:this.numOfCtScen
                ctScenId = this.ctScenProb(sCt,1);
                for sampleIx = 1:nScenarioSamples
                    rowIx = rowIx + 1;
                    linearMaskAll(rowIx,:) = [ctScenId, ...
                        scenarioIndexForSample(setupActive,sampleIx), ...
                        scenarioIndexForSample(rangeActive,sampleIx), ...
                        scenarioIndexForSample(gantryActive,sampleIx), ...
                        scenarioIndexForSample(couchActive,sampleIx)];
                end
            end

            %Mask for scenario selection
            if gantryActive || couchActive
                scenMaskSize = [this.numOfAvailableCtScen,this.totNumShiftScen, ...
                    this.totNumRangeScen,this.totNumGantryScen,this.totNumCouchScen];
                linearMask = linearMaskAll;
            else
                scenMaskSize = [this.numOfAvailableCtScen,this.totNumShiftScen, ...
                    this.totNumRangeScen];
                linearMask = linearMaskAll(:,1:3);
            end
            scenMask = false(scenMaskSize);
            maskSubscripts = mat2cell(linearMask, size(linearMask,1), ones(1,size(linearMask,2)));
            scenMask(sub2ind(scenMaskSize,maskSubscripts{:})) = true;
            totNumScen = sum(scenMask(:));

            scenarioValues = scenarios(:,2:end);
            scenProb = matRad_computeScenarioProbabilities(components,scenarioValues, ...
                this.ctScenProb,ctScenIds);
            [scenarioValues,ctScenIds,scenProb,scenWeight,scenarios,linearMask,scenMask] = ...
                matRad_filterZeroProbabilityScenarios(scenarioValues,ctScenIds, ...
                scenProb,scenWeight,scenarios,linearMask,scenMask);
            this.totNumScen = size(scenarioValues,1);
            totNumScen = sum(scenMask(:));
            this.setScenarioRealizations(components,scenarioValues,ctScenIds, ...
                scenProb,scenWeight,scenarios,linearMask,scenMask);

            if totNumScen ~= this.totNumScen
                matRad_cfg.dispWarning('Check Implementation of Total Scenario computation - given %d but found %d!',this.totNumScen,totNumScen);
                this.totNumScen = totNumScen;
            end
        end
    end


end

function n = scenarioCountForActiveDimension(isActive,nSamples)

if isActive
    n = nSamples;
else
    n = 1;
end

end

function scenIx = scenarioIndexForSample(isActive,sampleIx)

if isActive
    scenIx = sampleIx;
else
    scenIx = 1;
end

end
