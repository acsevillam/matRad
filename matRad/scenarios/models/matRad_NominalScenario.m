classdef matRad_NominalScenario < matRad_ScenarioModel
%  matRad_RandomScenarios
%  Implements a single nominal planning scenario
%
% constructor
%   matRad_NominalScenario()
%   matRad_NominalScenario(ct)
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
    
    properties (SetAccess = protected)
        name = 'nomScen';
    end

    methods
        function this = matRad_NominalScenario(ct)
            if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end            
            this@matRad_ScenarioModel(superclassArgs{:});

            this.initializeScenarioModel();
        end
        
        function scenarios = updateScenarios(this)
            this.numOfCtScen = size(this.ctScenProb,1);

            %set variables
            this.totNumShiftScen = 1;
            this.totNumRangeScen = 1;
            this.totNumGantryScen = 1;
            this.totNumCouchScen = 1;
            this.totNumScen = this.numOfCtScen; 

            components = this.getScenarioComponents();
            scenarioValues = zeros(this.numOfCtScen,numel(components));
            ctScenIds = this.ctScenProb(:,1);
            scenForProb = [ctScenIds scenarioValues];
            scenProb = matRad_computeScenarioProbabilities(components,scenarioValues, ...
                this.ctScenProb,ctScenIds);
            scenWeight = scenProb;

            %Mask for scenario selection
            scenMask = false(this.numOfAvailableCtScen,this.totNumShiftScen,this.totNumRangeScen);
            scenMask(ctScenIds,:,:) = true;

            %generic code
            [x{1}, x{2}, x{3}] = ind2sub(size(scenMask),find(scenMask));
            linearMask = cell2mat(x);
            [scenarioValues,ctScenIds,scenProb,scenWeight,scenForProb,linearMask,scenMask] = ...
                matRad_filterZeroProbabilityScenarios(scenarioValues,ctScenIds, ...
                scenProb,scenWeight,scenForProb,linearMask,scenMask);
            this.totNumScen = size(scenarioValues,1);
            totNumScen = sum(scenMask(:));
            this.setScenarioRealizations(components,scenarioValues,ctScenIds, ...
                scenProb,scenWeight,scenForProb,linearMask,scenMask);
            
            %Return variable
            scenarios = scenForProb;
            

            if totNumScen ~= this.totNumScen
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Check Implementation of Total Scenario computation - given %d but found %d!',this.totNumScen,totNumScen);
                this.totNumScen = totNumScen;
            end
        end

        function scenIx = sub2scenIx(this,ctScenInput,setupScenIx,rangeScenIx,ctScenReference)
            if nargin < 5
                scenIx = sub2scenIx@matRad_ScenarioModel(this,ctScenInput,setupScenIx,rangeScenIx);
            else
                scenIx = sub2scenIx@matRad_ScenarioModel(this,ctScenInput,setupScenIx,rangeScenIx,ctScenReference);
            end
        end
        
    end
end
