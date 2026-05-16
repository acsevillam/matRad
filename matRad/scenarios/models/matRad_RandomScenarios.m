classdef matRad_RandomScenarios < matRad_ScenarioModel
%  matRad_RandomScenarios
%  Implements randomly sampled scenarios
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
    
    properties
        includeNominalScenario = false; %Forces inclusion of the nominal scenario        
        nSamples = 10;                  %Standard number of random samples
    end

    %Deprecated Properties that were used
    properties (Dependent)
        numOfShiftScen;
        numOfRangeShiftScen;
    end

    properties (SetAccess=protected)
        name = 'Randomly sampled Scenarios';
        shortName   = 'rndScen';        
    end
    
    methods
        function this = matRad_RandomScenarios(ct)           
            if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end
            
            this@matRad_ScenarioModel(superclassArgs{:});

            %TODO: We could do this automatically in the superclass
            %Octave 5 has a bug there and throws an error
            this.updateScenarios();
        end

        %% Setters & Update
        function set.nSamples(this,nSamples)
            valid = isnumeric(nSamples) && isscalar(nSamples) && mod(nSamples,1) == 0 && nSamples > 0;
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for nSamples! Needs to be a positive integer!');
            end
            this.nSamples = nSamples;
            this.updateScenarios();
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

            % Multivariate normal sampling over active components only.
            scenarios = zeros(nScenarioSamples,numel(components));
            if any(activeIx)
                scales = [components(activeIx).scale];
                Sigma = diag(scales.^2);
                cs = chol(Sigma);
                scenarios(:,activeIx) = randn(nScenarioSamples,sum(activeIx)) * cs;
            end

            if this.includeNominalScenario
                % We include the nominal scenario by just replacing the
                % first one to keep the number of scenarios the same.
                scenarios(1,:) = 0;
            end

            % Handle 4D CT scenario ids.
            ctScenIds = repmat(this.ctScenProb(:,1)',size(scenarios,1),1);
            ctScenIds = ctScenIds(:);
            scenarios = horzcat(ctScenIds, repmat(scenarios,[this.numOfCtScen 1]));
            this.ctScenIx = ctScenIds;
            this.scenForProb = scenarios;

            this.scenProb = matRad_computeScenarioProbabilities(components,scenarios(:,2:end), ...
                this.ctScenProb,ctScenIds);

            % Scenario weight.
            this.scenWeight = [];
            for sCt = 1:this.numOfCtScen
                % Equal weights within a phase since they have been randomly
                % sampled; this is approximate when the nominal scenario is forced.
                this.scenWeight = [this.scenWeight; ...
                    this.ctScenProb(sCt,2) * ones(nScenarioSamples,1)./nScenarioSamples];
            end

            setupActive = any(strcmp({components.applicator},'setup') & [components.active]);
            rangeActive = any(strcmp({components.applicator},'range') & [components.active]);
            gantryActive = any(strcmp({components.applicator},'gantry') & [components.active]);
            couchActive = any(strcmp({components.applicator},'couch') & [components.active]);
            this.totNumShiftScen = scenarioCountForActiveDimension(setupActive,this.nSamples);
            this.totNumRangeScen = scenarioCountForActiveDimension(rangeActive,this.nSamples);
            this.totNumGantryScen = scenarioCountForActiveDimension(gantryActive,this.nSamples);
            this.totNumCouchScen = scenarioCountForActiveDimension(couchActive,this.nSamples);
            this.totNumScen = nScenarioSamples * this.numOfCtScen;

            this.relRangeShift = scenarios(:,6);
            this.absRangeShift = scenarios(:,5);
            this.isoShift = scenarios(:,2:4);
            this.maxAbsRangeShift = max(this.absRangeShift);
            this.maxRelRangeShift = max(this.relRangeShift);

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

            if gantryActive || couchActive
                scenMaskSize = [this.numOfAvailableCtScen,this.totNumShiftScen, ...
                    this.totNumRangeScen,this.totNumGantryScen,this.totNumCouchScen];
                this.linearMask = linearMaskAll;
            else
                scenMaskSize = [this.numOfAvailableCtScen,this.totNumShiftScen, ...
                    this.totNumRangeScen];
                this.linearMask = linearMaskAll(:,1:3);
            end

            this.scenMask = false(scenMaskSize);
            maskSubscripts = mat2cell(this.linearMask, size(this.linearMask,1), ...
                ones(1,size(this.linearMask,2)));
            this.scenMask(sub2ind(scenMaskSize,maskSubscripts{:})) = true;
            totNumScen = sum(this.scenMask(:));
            this.finalizeScenarioRealizations();

            if totNumScen ~= this.totNumScen
                matRad_cfg.dispWarning('Check Implementation of Total Scenario computation - given %d but found %d!',this.totNumScen,totNumScen);
                this.totNumScen = totNumScen;
            end
        end
    end

    methods (Access = protected)

        function validateScenarioDimensionSupport(~,~)
            % Random sampling supports all currently declared continuous
            % scenario dimensions.
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
