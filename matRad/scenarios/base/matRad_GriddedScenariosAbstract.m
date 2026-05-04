classdef (Abstract) matRad_GriddedScenariosAbstract < matRad_ScenarioModel
% matRad_GriddedScenariosAbstract base class for grid-based scenario models
%
% This abstract class provides common setup and range grid construction for
% worst-case, importance, and truncated-importance scenario models. Setup
% shifts are represented in mm, absolute range shifts in mm, and relative
% range shifts as fractions internally after conversion from percent input.
%
% Usage:
%   this = matRad_GriddedScenariosAbstract()
%   this = matRad_GriddedScenariosAbstract(ct)
%
% input
%   ct:                 matRad ct struct used to derive available CT
%                       scenario ids and CT scenario probabilities
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
        %includeNominalScenario = true;        
        combinations = 'none'; % Can be 'none', 'shift', or 'all' to control scenario combinations
        combineRange = true; % Whether to treat absolute and relative range as one shift or separate scenarios
    end
    
    %Each subclass needs to define how many gridpoints it uses and if this
    %can be set or not
    properties (Abstract)
        numOfSetupGridPoints;
        numOfRangeGridPoints;
    end

    properties (Constant)
        validCombinationTypes = {'all','none','shift'};
    end

    methods
        function this = matRad_GriddedScenariosAbstract(ct)      
                         
            if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end
            
            this@matRad_ScenarioModel(superclassArgs{:});
        end
        
        function set.combineRange(this,combineRange_)
            valid = isscalar(combineRange_) && (isnumeric(combineRange_) || islogical(combineRange_));
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for combineRange! Needs to be a boolean / logical value!');
            end
            this.combineRange = combineRange_;
            this.requestScenarioUpdate();
        end

        %% set methods
        %{
        function set.includeNominalScenario(this,includeNomScen)
            valid = isscalar(includeNomScen) && (isnumeric(includeNomScen) || islogical(includeNomScen));
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for includeNominalScenario! Needs to be a boolean / logical value!');
            end
            this.includeNominalScenario = includeNomScen;
            this.updateScenarios();
        end
        %}
    
        function set.combinations(this,combinations_)
            valid = any(strcmp(combinations_,this.validCombinationTypes));
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for combinations! Needs to be one of the strings %s!',strjoin(this.validCombinationTypes,' / '));
            end
            this.combinations = combinations_;
            this.requestScenarioUpdate();
        end

        function scenarios = updateScenarios(this)
            matRad_cfg = MatRad_Config.instance();

            %
            this.numOfCtScen = size(this.ctScenProb,1);
            components = this.getScenarioComponents();
            componentScales = [components.scale];
            componentScales(~[components.active]) = 0;

            %Get the maximum, i.e., worst case shifts
            wcSetupShifts = this.wcSigma * componentScales(1:3);
            
            %% Create gridded setup shifts
            %Create grid vectors for setup shifts
            setupShiftGrid = zeros(this.numOfSetupGridPoints,numel(wcSetupShifts));
            %{
            if mod(this.numOfSetupGridPoints,2) == 0 && this.includeNominalScenario
                matRad_cfg.dispWarning('Obtaining Setup Shifts: Including the nominal scenario with even number of grid points creates asymmetrical shifts!');
            end
            %}

            for i = 1:numel(wcSetupShifts)
                setupShiftGrid(:,i) = linspace(-wcSetupShifts(i),wcSetupShifts(i),this.numOfSetupGridPoints);
                %{
                if this.includeNominalScenario 
                      
                    [~,ix] = min(abs(setupShiftGrid(:,i)));
                    setupShiftGrid(ix,i) = 0;    
                end  
                %}
            end
            
            %Now create vector of all shifts for different combinatorial
            %settings
            switch this.combinations
                case 'none'
                    %Create independent shifts
                    griddedSetupShifts = [];
                    for i=1:size(setupShiftGrid,2)
                        tmpGrid = zeros(size(setupShiftGrid,1),3);
                        tmpGrid(:,i) = setupShiftGrid(:,i);
                        griddedSetupShifts = [griddedSetupShifts; tmpGrid];
                    end                    
                case {'shift','all'}
                    [X,Y,Z] = meshgrid(setupShiftGrid(:,1),setupShiftGrid(:,2),setupShiftGrid(:,3));
                    griddedSetupShifts = [X(:), Y(:), Z(:)];    
                otherwise
                    matRad_cfg.dispError('Invalid value for combinations! This sanity check should never be reached!');
            end

            griddedSetupShifts = matRad_ImportanceScenarios.uniqueStableRowsCompat(griddedSetupShifts);
            shiftNomScenIx = find(all(griddedSetupShifts == zeros(1,3),2));            
            
            if ~isempty(shiftNomScenIx) %|| this.includeNominalScenario
                if ~isempty(shiftNomScenIx)
                    griddedSetupShifts(shiftNomScenIx,:) = [];
                end
                griddedSetupShifts = [0 0 0; griddedSetupShifts];
            end
                        
            this.totNumShiftScen = size(griddedSetupShifts,1);
                                
            %% Create gridded range shifts
            %Obtain worst case range shifts
            wcRangeShifts = this.wcSigma * componentScales(4:5);
            
            rangeShiftGrid = zeros(this.numOfRangeGridPoints,numel(wcRangeShifts));  
            %{
            if mod(this.numOfRangeGridPoints,2) == 0 && this.includeNominalScenario
                matRad_cfg.dispWarning('Obtaining Range Shifts: Including the nominal scenario with even number of grid points creates asymmetrical shifts!');
            end
            %}

            for i = 1:numel(wcRangeShifts)
                rangeShiftGrid(:,i) = linspace(-wcRangeShifts(i),wcRangeShifts(i),this.numOfRangeGridPoints);
                
                %{
                if this.includeNominalScenario 
                    [~,ix] = min(abs(rangeShiftGrid(:,i)));
                    rangeShiftGrid(ix,i) = 0;    
                end
                %}
            end

            if this.combineRange
                griddedRangeShifts = rangeShiftGrid;
            else                
                [rngAbs,rngRel] = meshgrid(rangeShiftGrid(:,1),rangeShiftGrid(:,2));
                griddedRangeShifts = [rngAbs(:),rngRel(:)];
            end

            %Remove duplicate scenarios and update number of shifts
            griddedRangeShifts = this.uniqueStableRowsCompat(griddedRangeShifts); 

            rangeNomScenIx = find(all(griddedRangeShifts == zeros(1,2),2));            
            
            if ~isempty(rangeNomScenIx) %|| this.includeNominalScenario
                if ~isempty(rangeNomScenIx)
                    griddedRangeShifts(rangeNomScenIx,:) = [];
                end
                griddedRangeShifts = [0 0; griddedRangeShifts];
            end

            this.totNumRangeScen = size(griddedRangeShifts,1);
            this.totNumGantryScen = 1;
            this.totNumCouchScen = 1;
                       
            %Aggregate scenarios
            switch this.combinations
                case {'none','shift'}
                    scenarios = zeros(this.totNumShiftScen + this.totNumRangeScen,5);
                    scenarios(1:this.totNumShiftScen,1:3) = griddedSetupShifts;
                    scenarios(this.totNumShiftScen+1:this.totNumShiftScen+this.totNumRangeScen,4:5) = griddedRangeShifts;

                    %create the linear mask of scenarios
                    linearMaskTmp = ones(size(scenarios,1),3);
                    linearMaskTmp(1:this.totNumShiftScen,2) = (1:this.totNumShiftScen)';
                    linearMaskTmp(this.totNumShiftScen+1:this.totNumShiftScen+this.totNumRangeScen,3) = (1:this.totNumRangeScen)';

                    [scenarios,ia] = matRad_ImportanceScenarios.uniqueStableRowsCompat(scenarios);
                    linearMaskTmp = linearMaskTmp(ia,:);

                case 'all'
                    %Prepare scenario matrix by replicating shifts
                    %with the number of range scenarios
                    scenarios = repmat(griddedSetupShifts,this.totNumRangeScen,1);
                    scenarios = [scenarios zeros(size(scenarios,1),2)];
                    
                    %create the linear mask of scenarios
                    linearMaskTmp = ones(size(scenarios,1),3);
                    for r = 1:this.totNumRangeScen                    
                        offset = (r-1)*this.totNumShiftScen;
                        ixR = (offset + 1):(offset + this.totNumShiftScen);
                        scenarios(ixR,4:5) = repmat(griddedRangeShifts(r,:),this.totNumShiftScen,1);

                        %Set linear mask
                        linearMaskTmp(ixR,2) = (1:this.totNumShiftScen)';
                        linearMaskTmp(ixR,3) = r;
                    end


                    %create the linear mask of scenarios
                    [scenarios,ia] = matRad_ImportanceScenarios.uniqueStableRowsCompat(scenarios);
                    linearMaskTmp = linearMaskTmp(ia,:);
                otherwise
                    matRad_cfg.dispError('Invalid value for combinations! This sanity check should never be reached!');
            end

            %if ~this.includeNominalScenario
            %    nomScen = all(scenarios == zeros(1,5),2);
            %    scenarios(nomScen,:) = [];
            %    linearMaskTmp(nomScen,:) = [];
            %end

            %Handle 4D CT scenario ids
            ctScenIds = repmat(this.ctScenProb(:,1)',size(scenarios,1),1);
            ctScenIds = ctScenIds(:);
            scenarios = horzcat(ctScenIds, repmat(scenarios,[this.numOfCtScen 1]));
            linearMaskTmp = repmat(linearMaskTmp,this.numOfCtScen,1);
            linearMaskTmp(:,1) = ctScenIds;
            %Finalize meta information
            this.totNumScen = size(scenarios,1);

            scenMask = false(this.numOfAvailableCtScen,this.totNumShiftScen,this.totNumRangeScen);
            
            maskIx = sub2ind(size(scenMask),linearMaskTmp(:,1),linearMaskTmp(:,2),linearMaskTmp(:,3));
            scenMask(maskIx) = true;

            scenarioValues = scenarios(:,2:end);
            scenProb = matRad_computeScenarioProbabilities(components,scenarioValues, ...
                this.ctScenProb,ctScenIds);
            scenWeight = scenProb;
            [scenarioValues,ctScenIds,scenProb,scenWeight,scenarios,linearMaskTmp,scenMask] = ...
                matRad_filterZeroProbabilityScenarios(scenarioValues,ctScenIds, ...
                scenProb,scenWeight,scenarios,linearMaskTmp,scenMask);
            this.totNumScen = size(scenarioValues,1);
            this.setScenarioRealizations(components,scenarioValues,ctScenIds, ...
                scenProb,scenWeight,scenarios,linearMaskTmp,scenMask);
        end
    end

    methods (Access = protected)
        function validateScenarioDimensionSupport(~,scenarioDimensionActive)
            if any(strcmp(scenarioDimensionActive,'gantry')) || any(strcmp(scenarioDimensionActive,'couch'))
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError(['Gantry/couch uncertainty dimensions are only supported by random scenario models for now. ' ...
                    'Use matRad_RandomScenarios (rndScen) for angular uncertainty.']);
            end
        end
    end

    methods (Static)
        function [uniqueStableRows,ia] = uniqueStableRowsCompat(values)
            %This is a compatability wrapper to call unique without sorting            
            
            matRad_cfg = MatRad_Config.instance();
            
            if matRad_cfg.isOctave
                %https://stackoverflow.com/questions/37828894/
                [~,ia,~] = unique(values,'rows','first');
                ia = sort(ia);
                uniqueStableRows = values(ia,:);
            else
                [uniqueStableRows,ia] = unique(values,'rows','stable');
            end
        end
    end
end
