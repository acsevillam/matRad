classdef (Abstract) matRad_ScenarioModel < handle
%  matRad_ScenarioModel
%  This is an abstract interface class to define Scenario Models for use in
%  robust treatment planning and uncertainty analysis.
%  Subclasses should at least implement the update() function to generate
%  their own scenarios.
%
% constructor (Abstract)
%   matRad_ScenarioModel()
%   matRad_ScenarioModel(ct)
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
    
    properties (AbortSet = true) %We use AbortSet = true here to avoid updates when 
        %Uncertainty model
        rangeRelSD  = 3.5;                % given in %
        rangeAbsSD  = 1;                  % given in [mm]
        shiftSD     = [2.25 2.25 2.25];   % given in [mm]
        gantryAngleSD = 1;                % given in [deg]
        couchAngleSD  = 1;                % given in [deg]
        numOfBeams = 0;                   % number of beams for per-beam angular uncertainty components
        wcSigma     = 1;                  % Multiplier to compute the worst case / maximum shifts
        scenarioDimensionActive = {'ct','setup','range'}; % active uncertainty dimensions/applicators

        ctScenProb  = [1 1];              % Ct Scenarios to be included in the model. Left column: Scenario Index. Right column: Scenario Probability        
    end

    properties (Abstract,SetAccess = protected)
        name
        shortName
    end

    properties (Dependent)
        wcFactor;
    end
   
    properties (SetAccess = protected)
        numOfCtScen;            % total number of CT scenarios used
        numOfAvailableCtScen;   % total number of CT scenarios existing in ct structure
        ctScenIx;               % map of all ct scenario indices per scenario


        % these parameters will be filled according to the choosen scenario type
        isoShift;
        relRangeShift;
        absRangeShift;
        gantryAngleOffset;
        couchAngleOffset;

        maxAbsRangeShift;
        maxRelRangeShift;
        
        totNumShiftScen;        % total number of shift scenarios in x,y and z direction
        totNumRangeScen;        % total number of range and absolute range scenarios
        totNumGantryScen;       % total number of gantry angle scenarios
        totNumCouchScen;        % total number of couch angle scenarios
        totNumScen;             % total number of samples 
        
        scenForProb;            % matrix for probability calculation - each row denotes one scenario, whereas columns denotes the realization value
        scenProb;               % probability of each scenario stored in a vector (according to uncertainty model)
        scenWeight;             % weight of scenario relative to the underlying uncertainty model (depends on how scenarios are chosen / sampled)
        scenMask;
        linearMask;

        scenarioComponents;     % scenario realization components grouped by applicator
        scenarioValueNames;     % ordered names for columns in scenarioValues
        scenarioValues;         % scenario realizations, one row per scenario
        scenarioIdList;         % stable public scenario ids
        scenarioCtScenIds;      % ct scenario id for each realization row
        scenarioApplicators;    % applicator objects supported by this model
    end
    
    methods
        function this = matRad_ScenarioModel(ct)
            if nargin == 0 || isempty(ct)
                this.numOfCtScen = 1;
                this.numOfAvailableCtScen = 1;
            else
                this.numOfCtScen = ct.numOfCtScen;
                this.numOfAvailableCtScen = ct.numOfCtScen;
            end

            this.ctScenProb = [(1:this.numOfCtScen)', ones(this.numOfCtScen,1)./this.numOfCtScen]; %Equal probability to be in each phase of the 4D ct
            
            %TODO: We could do this here automatically in the constructor, but
            %Octave 5 has a bug here and throws an error
            %this.updateScenarios();
        end

        function listAllScenarios(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('Listing all scenarios...\n');
            matRad_cfg.dispInfo('\t#\tctScen\txShift\tyShift\tzShift\tabsRng\trelRng\tprob.\n');
            for s = 1:size(this.scenForProb,1)
                str = [num2str(this.scenForProb(s,1),'%d\t'),sprintf('\t'), num2str(this.scenForProb(s,2:end),'\t%.3f')];
                matRad_cfg.dispInfo('\t%d\t%s\t%.3f\n',s,str,this.scenProb(s));
            end
        end

        %% SETTERS & UPDATE
        function set.rangeRelSD(this,rangeRelSD)
            valid = isnumeric(rangeRelSD) && isscalar(rangeRelSD) && rangeRelSD >= 0;
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for rangeRelSD! Needs to be a real positive scalar!');
            end
            this.rangeRelSD = rangeRelSD;
            this.updateScenarios();
        end

        function set.rangeAbsSD(this,rangeAbsSD)
            valid = isnumeric(rangeAbsSD) && isscalar(rangeAbsSD) && rangeAbsSD >= 0;
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for rangeAbsSD! Needs to be a real positive scalar!');
            end
            this.rangeAbsSD = rangeAbsSD;
            this.updateScenarios();
        end

        function set.shiftSD(this,shiftSD)
            valid = isnumeric(shiftSD) && isrow(shiftSD) && numel(shiftSD) == 3 && all(shiftSD > 0);
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for shiftSD! Needs to be 3-element numeric row vector!');
            end
            this.shiftSD = shiftSD;
            this.updateScenarios();
        end

        function set.gantryAngleSD(this,gantryAngleSD)
            valid = isnumeric(gantryAngleSD) && isscalar(gantryAngleSD) && gantryAngleSD >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for gantryAngleSD! Needs to be a real non-negative scalar!');
            end
            this.gantryAngleSD = gantryAngleSD;
            this.updateScenarios();
        end

        function set.couchAngleSD(this,couchAngleSD)
            valid = isnumeric(couchAngleSD) && isscalar(couchAngleSD) && couchAngleSD >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for couchAngleSD! Needs to be a real non-negative scalar!');
            end
            this.couchAngleSD = couchAngleSD;
            this.updateScenarios();
        end

        function set.numOfBeams(this,numOfBeams)
            valid = isnumeric(numOfBeams) && isscalar(numOfBeams) && isfinite(numOfBeams) && ...
                round(numOfBeams) == numOfBeams && numOfBeams >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for numOfBeams! Needs to be a non-negative integer scalar!');
            end
            this.numOfBeams = numOfBeams;
            this.updateScenarios();
        end

        function set.wcSigma(this,wcSigma)
            valid = isnumeric(wcSigma) && isscalar(wcSigma) && wcSigma >= 0;
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for wcSigma! Needs to be a real positive scalar!');
            end
            this.wcSigma = wcSigma;
            this.updateScenarios();
        end

        function set.ctScenProb(this,ctScenProb)
            valid = isnumeric(ctScenProb) && ismatrix(ctScenProb) && size(ctScenProb,2) == 2 && all(round(ctScenProb(:,1)) == ctScenProb(:,1)) && all(ctScenProb(:) >= 0);
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for used ctScenProb! Needs to be a valid 2-column matrix with left column representing the scenario index and right column representing the appropriate probabilities [0,1]!');
            end            
            this.ctScenProb = ctScenProb;
            this.updateScenarios();
        end

        function set.scenarioDimensionActive(this,scenarioDimensionActive)
            scenarioDimensionActive = matRad_normalizeScenarioDimensionActive(scenarioDimensionActive);
            this.validateScenarioDimensionSupport(scenarioDimensionActive);
            this.scenarioDimensionActive = scenarioDimensionActive;
            this.updateScenarios();
        end


        function scenarios = updateScenarios(this)            
            %This function will always update the scenarios given the
            %current property settings

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('This abstract function needs to be implemented!');
        end

        function newInstance = extractSingleScenario(this,scenNum)
            newInstance = matRad_NominalScenario();
            
            ctScenNum = this.linearMask(scenNum,1);
            
            %First set properties that force an update
            newInstance.numOfCtScen         = 1;            
            newInstance.ctScenProb          = this.ctScenProb(ctScenNum,:);

            %Now overwrite existing variables for correct probabilties and
            %error realizations
            newInstance.scenForProb         = this.scenForProb(scenNum,:);
            newInstance.relRangeShift       = this.scenForProb(scenNum,6);
            newInstance.absRangeShift       = this.scenForProb(scenNum,5);
            newInstance.isoShift            = this.scenForProb(scenNum,2:4);
            newInstance.scenProb            = this.scenProb(scenNum);
            newInstance.scenWeight          = this.scenWeight(scenNum);
            newInstance.maxAbsRangeShift    = max(abs(this.absRangeShift(scenNum)));
            newInstance.maxRelRangeShift    = max(abs(this.relRangeShift(scenNum)));
            newInstance.scenMask            = false(this.numOfAvailableCtScen,1,1);
            newInstance.linearMask          = [newInstance.ctScenIx 1 1];
            
            newInstance.scenMask(newInstance.linearMask(:,1),newInstance.linearMask(:,2),newInstance.linearMask(:,3)) = true;
            newInstance.finalizeScenarioRealizations();
            %newInstance.updateScenarios();
        end
        
        function scenIx = sub2scenIx(this,ctScen,shiftScen,rangeShiftScen)
            %Returns linear index in the scenario cell array from scenario
            %subscript indices
            if ~isvector(this.scenMask)
                scenIx = sub2ind(size(this.scenMask),ctScen,shiftScen,rangeShiftScen);
            else
                scenIx = this.ctScenIx(ctScen);
            end
        end

        function scenNum = scenNum(this,fullScenIx)
            %gets number of scneario from full scenario index in scenMask
            scenNum = find(find(this.scenMask) == fullScenIx);
        end

        function ids = scenarioIds(this)
            ids = this.scenarioIdList(:);
        end

        function n = numScenarios(this)
            n = numel(this.scenarioIdList);
        end

        function scenario = getScenario(this,scenarioId)
            scenarioRowIx = this.resolveScenarioRowIx(scenarioId);
            scenario = struct();
            scenario.id = this.scenarioIdList(scenarioRowIx);
            scenario.ctScenId = this.scenarioCtScenIds(scenarioRowIx);
            scenario.values = this.rowToValueStruct(this.scenarioValues(scenarioRowIx,:));
            scenario.probability = this.scenProb(scenarioRowIx);
            scenario.weight = this.scenWeight(scenarioRowIx);
        end

        function ctScenId = getCtScenario(this,scenarioId)
            scenarioRowIx = this.resolveScenarioRowIx(scenarioId);
            ctScenId = this.scenarioCtScenIds(scenarioRowIx);
        end

        function value = getValue(this,scenarioId,componentName)
            scenarioRowIx = this.resolveScenarioRowIx(scenarioId);
            componentIx = this.findScenarioComponentIndex(componentName);
            value = this.scenarioValues(scenarioRowIx,componentIx);
        end

        function values = getValues(this,scenarioId,componentNames)
            if ischar(componentNames)
                componentNames = {componentNames};
            end

            values = zeros(1,numel(componentNames));
            for i = 1:numel(componentNames)
                values(i) = this.getValue(scenarioId,componentNames{i});
            end
        end

        function shift = getSetupShift(this,scenarioId)
            shift = this.getValues(scenarioId,{'setup.x','setup.y','setup.z'});
        end

        function rangeShift = getRangeShift(this,scenarioId)
            rangeShift = this.getValues(scenarioId,{'range.absolute','range.relative'});
        end

        function gantryOffsets = getGantryAngleOffset(this,scenarioId)
            gantryOffsets = this.getBeamAngleOffsets(scenarioId,'gantry');
        end

        function couchOffsets = getCouchAngleOffset(this,scenarioId)
            couchOffsets = this.getBeamAngleOffsets(scenarioId,'couch');
        end

        function tf = hasActiveScenarioDimension(this,dimensionName)
            if isstring(dimensionName) && isscalar(dimensionName)
                dimensionName = char(dimensionName);
            end
            tf = any(strcmp(this.scenarioDimensionActive,dimensionName));
        end

        function tf = hasActiveAngularScenarioDimension(this)
            tf = this.hasActiveScenarioDimension('gantry') || ...
                this.hasActiveScenarioDimension('couch');
        end

        function stf = applyScenarioToStf(this,scenarioId,stf)
            applicators = this.getScenarioApplicators();
            for i = 1:numel(applicators)
                applicator = applicators{i};
                if ismethod(applicator,'applyToStf')
                    stf = applicator.applyToStf(this,scenarioId,stf);
                end
            end
        end

        function ids = getNominalScenarioIds(this)
            ids = this.scenarioIdList(all(abs(this.scenarioValues) <= eps,2));
        end

        function fullScenIx = getDijScenarioIndex(this,scenarioId)
            scenarioRowIx = this.resolveScenarioRowIx(scenarioId);
            fullScenIxs = find(this.scenMask);
            fullScenIx = fullScenIxs(scenarioRowIx);
        end

        function fullScenIx = getDijScenarioIndexBySubscripts(this,ctScenInput,shiftScenIx,rangeScenIx,ctScenReference)
            if nargin < 5 || isempty(ctScenReference)
                ctScenReference = 'position';
            end

            ctScenId = helper_resolveCtScenarioId(this,ctScenInput,ctScenReference);
            helper_validatePositiveIntegerScalar(shiftScenIx,'shiftScenIx');
            helper_validatePositiveIntegerScalar(rangeScenIx,'rangeScenIx');

            scenMaskSize = size(this.scenMask);
            scenMaskSize(end + 1:3) = 1;
            if ctScenId > scenMaskSize(1) || shiftScenIx > scenMaskSize(2) || rangeScenIx > scenMaskSize(3)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Scenario subscript exceeds the scenario mask dimensions.');
            end
            fullScenIx = sub2ind(scenMaskSize,ctScenId,shiftScenIx,rangeScenIx);
        end

        function tf = isScenarioActiveBySubscripts(this,ctScenInput,shiftScenIx,rangeScenIx,ctScenReference)
            if nargin < 5
                ctScenReference = 'position';
            end
            fullScenIx = this.getDijScenarioIndexBySubscripts(ctScenInput,shiftScenIx,rangeScenIx,ctScenReference);
            tf = this.scenMask(fullScenIx);
        end

        function scenarioRowIx = getScenarioRowIndexFromDijIndex(this,fullScenIx)
            scenarioRowIx = find(find(this.scenMask) == fullScenIx,1,'first');
        end

        function applicators = getScenarioApplicators(this)
            applicators = this.scenarioApplicators;
        end

        function mask = getDijActiveMask(this)
            mask = this.scenMask;
        end

        function sz = getDijContainerSize(this)
            sz = size(this.scenMask);
        end
        
        %% Deprecated functions / properties
        function newInstance = extractSingleNomScen(this,~,scenIdx)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The function extractSingleNomScen of the scenario class will soon be deprecated! Use extractSingleScenario instead!');
            newInstance = this.extractSingleScenario(scenIdx);
        end

        function t = TYPE(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property TYPE of the scenario class will soon be deprecated!');
            t = this.shortName;
        end

        function value = get.wcFactor(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property wcFactor of the scenario class will soon be deprecated!');
            value = this.wcSigma;
        end

        function set.wcFactor(this,value)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property wcFactor of the scenario class will soon be deprecated!');
            this.wcSigma = value;
        end

    end

    methods (Access = protected)
        function finalizeScenarioRealizations(this)
            components = this.getScenarioComponents();
            scenarioValues = this.scenForProb(:,2:end);
            if numel(components) ~= size(scenarioValues,2)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError(['Scenario realization values do not match the declared uncertainty ' ...
                    'components. Found %d value columns for %d components.'], ...
                    size(scenarioValues,2), numel(components));
            end

            this.scenarioComponents = components;
            this.scenarioValueNames = {components.name};
            this.scenarioValues = scenarioValues;
            this.scenarioIdList = (1:size(scenarioValues,1))';
            this.scenarioCtScenIds = this.scenForProb(:,1);
            this.scenarioApplicators = { ...
                matRad_CtScenarioApplicator(), ...
                matRad_SetupShiftApplicator(), ...
                matRad_RangeShiftApplicator(), ...
                matRad_GantryAngleApplicator(), ...
                matRad_CouchAngleApplicator()};

            this.gantryAngleOffset = this.extractScenarioColumns(this.beamAngleComponentNames('gantry'));
            this.couchAngleOffset = this.extractScenarioColumns(this.beamAngleComponentNames('couch'));
        end

        function components = getScenarioComponents(this)
            components = matRad_createScenarioComponents(this.shiftSD,this.rangeAbsSD, ...
                this.rangeRelSD,this.scenarioDimensionActive,this.numOfBeams, ...
                this.gantryAngleSD,this.couchAngleSD);
        end

        function validateScenarioDimensionSupport(~,scenarioDimensionActive)
            if any(strcmp(scenarioDimensionActive,'gantry')) || any(strcmp(scenarioDimensionActive,'couch'))
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError(['Gantry/couch uncertainty dimensions are declared in the scenario architecture ' ...
                    'but are not supported by this scenario model. Use an angular-capable scenario model instead.']);
            end
        end

        function componentNames = beamAngleComponentNames(this,applicatorName)
            componentNames = cell(1,this.numOfBeams);
            for i = 1:this.numOfBeams
                componentNames{i} = sprintf('%s.beam%d',applicatorName,i);
            end
        end

        function offsets = getBeamAngleOffsets(this,scenarioId,applicatorName)
            componentNames = this.beamAngleComponentNames(applicatorName);
            if isempty(componentNames) || ~this.hasActiveScenarioDimension(applicatorName)
                offsets = zeros(1,this.numOfBeams);
            else
                offsets = this.getValues(scenarioId,componentNames);
            end
        end

        function scenarioRowIx = resolveScenarioRowIx(this,scenarioId)
            scenarioRowIx = find(this.scenarioIdList == scenarioId,1,'first');
            if isempty(scenarioRowIx)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Scenario id %d does not exist in this scenario model.',scenarioId);
            end
        end

        function componentIx = findScenarioComponentIndex(this,componentName)
            if isstring(componentName) && isscalar(componentName)
                componentName = char(componentName);
            end

            componentIx = find(strcmp(this.scenarioValueNames,componentName),1,'first');
            if isempty(componentIx)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Scenario component "%s" does not exist in this scenario model.',componentName);
            end
        end

        function valueStruct = rowToValueStruct(this,valueRow)
            valueStruct = struct();
            for i = 1:numel(this.scenarioValueNames)
                fieldName = matlab.lang.makeValidName(this.scenarioValueNames{i});
                valueStruct.(fieldName) = valueRow(i);
            end
        end

        function values = extractScenarioColumns(this,componentNames)
            values = zeros(size(this.scenarioValues,1),numel(componentNames));
            for i = 1:numel(componentNames)
                componentIx = find(strcmp(this.scenarioValueNames,componentNames{i}),1,'first');
                if ~isempty(componentIx)
                    values(:,i) = this.scenarioValues(:,componentIx);
                end
            end
        end
    end

    methods (Static)
 
        function classList = getAvailableModels()
            matRad_cfg = MatRad_Config.instance();
            
            % Use the scenarios root folder so models can live in a
            % dedicated subfolder while keeping the public class names.
            scenariosRoot = fileparts(fileparts(mfilename('fullpath')));
            folders = {scenariosRoot};
            folders = [folders matRad_cfg.userfolders];

            persistent metaScenarioModels lastOptionalPaths
            
            %First we do a sanity check if persistently stored metaclasses are valid
            if ~matRad_cfg.isOctave && ~isempty(metaScenarioModels) && ~all(cellfun(@isvalid,metaScenarioModels))
                matRad_cfg.dispWarning('Found invalid ScenarioModels, updating model cache.');
                metaScenarioModels = [];
            end

            if isempty(metaScenarioModels) || (~isempty(lastOptionalPaths) && ~isequal(lastOptionalPaths, folders))
                lastOptionalPaths = folders;
                metaScenarioModels = matRad_findSubclasses(meta.class.fromName(mfilename('class')),'folders',folders,'includeSubfolders',true);
            end
            classList = matRad_identifyClassesByConstantProperties(metaScenarioModels,'shortName','defaults',{'nomScen'});

            if isempty(classList)
                matRad_cfg.dispError('No models found in paths %s',strjoin(folders,'\n'));
            end
        end
        
        function model = create(modelMetadata,ct)
            if isa(modelMetadata,'matRad_ScenarioModel')
                model = modelMetadata;
                return;
            end
            
            matRad_cfg = MatRad_Config.instance();

            if ischar(modelMetadata) || isstring(modelMetadata)
                modelMetadata = struct('model',modelMetadata);
            end

            modelClassList = matRad_ScenarioModel.getAvailableModels();
            modelNames = {modelClassList.shortName};
            
            if ~isfield(modelMetadata,'model') || ~any(strcmp(modelNames,modelMetadata.model))
                matRad_cfg.dispWarning('Scenario Model not found, creating nominal scenario instead!');
                modelMetadata.model = 'nomScen';
            end

            usedModel = find(strcmp(modelNames,modelMetadata.model));
            
            if ~isscalar(usedModel)
                usedModel = usedModel(1);
            end
            
            modelClassInfo = modelClassList(usedModel);

            if nargin < 2
                model = modelClassInfo.handle();
            else
                model = modelClassInfo.handle(ct);
            end

            modelMetadata = rmfield(modelMetadata,'model');
            
            %Now overwrite properties
            fields = fieldnames(modelMetadata);
            
            % iterate over all fieldnames and try to set the
            % corresponding properties inside the engine
            if matRad_cfg.isOctave
                c2sWarningState = warning('off','Octave:classdef-to-struct');                
            end
            
            for i = 1:length(fields)
                try
                    field = fields{i};
                    if helper_isDerivedScenarioMetadataField(field)
                        continue
                    end

                    if matRad_ispropCompat(model,field)
                        model.(field) = matRad_recursiveFieldAssignment(model.(field),modelMetadata.(field),true);
                    else
                        matRad_cfg.dispWarning('Not able to assign property ''%s'' from scenario model struct!',field);
                    end
                catch ME
                    % catch exceptions when the model has no properties,
                    % which are defined in the struct.
                    % When defining an engine with custom setter and getter
                    % methods, custom exceptions can be caught here. Be
                    % careful with Octave exceptions!
                    matRad_cfg = MatRad_Config.instance();
                    switch ME.identifier
                        case 'MATLAB:noPublicFieldForClass'
                            matRad_cfg.dispWarning('Not able to assign property from scenario model struct: %s',ME.message);
                        otherwise
                            matRad_cfg.dispWarning('Problem while setting up scenario Model from struct:%s %s',field,ME.message);
                    end
                end
            end
            
            if matRad_cfg.isOctave
                warning(c2sWarningState.state,'Octave:classdef-to-struct');                
            end
        end
    end
end

function tf = helper_isDerivedScenarioMetadataField(field)

tf = any(strcmp(field,{'scenarioComponents','scenarioValueNames','scenarioValues', ...
                       'scenarioIdList','scenarioCtScenIds','scenarioApplicators'}));

end

function ctScenId = helper_resolveCtScenarioId(scenarioModel,ctScenInput,ctScenReference)

helper_validatePositiveIntegerScalar(ctScenInput,'ctScen');

if isstring(ctScenReference) && isscalar(ctScenReference)
    ctScenReference = char(ctScenReference);
end

if ~ischar(ctScenReference)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('ctScenReference must be ''position'' or ''id''.');
end

switch lower(ctScenReference)
    case {'position','ctscenposition'}
        if ctScenInput > size(scenarioModel.ctScenProb,1)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('CT scenario position %d exceeds the scenario model size.',ctScenInput);
        end
        ctScenId = scenarioModel.ctScenProb(ctScenInput,1);
    case {'id','ctscenid'}
        ctScenId = ctScenInput;
        if ~any(scenarioModel.ctScenProb(:,1) == ctScenId)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Could not find CT scenario %d in the scenario model.',ctScenId);
        end
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('ctScenReference must be ''position'' or ''id''.');
end

end

function helper_validatePositiveIntegerScalar(value,valueName)

if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        round(value) == value && value >= 1)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('%s must be a positive integer scalar.',valueName);
end

end
