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
% input
%   ct:                 ct cube
%
% References
%   -
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
        ctScenProb  = [1 1];              % CT scenarios included in the model. Left column: CT scenario id. Right column: CT scenario probability
    end

    properties (Abstract,SetAccess=protected)
        name
    end

    properties (Dependent, Hidden)
        wcFactor;
        ctScenIx;               % deprecated: use scenarioCtScenIds
        scenarioCtScen;         % deprecated: use scenarioCtScenIds
    end

    properties (SetAccess = protected)
        numOfCtScen;            % total number of CT scenarios used
        numOfAvailableCtScen;   % total number of CT scenarios existing in ct structure


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

    properties (Access = private)
        updateScenariosOnChange = false;
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
        end

        function listAllScenarios(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('Listing all scenarios...\n');
            matRad_cfg.dispInfo('\t#\txShift\tyShift\tzShift\tabsRng\trelRng\tprob.\n');
            for s = 1:this.numScenarios()
                str = num2str([this.getCtScenario(s) this.getSetupShift(s) this.getRangeShift(s)],'\t%.3f');
                matRad_cfg.dispInfo('\t%d\t%s\t%.3f\n',s,str,this.scenProb(s));
            end
        end

        %% SETTERS & UPDATE
        function set.rangeRelSD(this,rangeRelSD)
            valid = isnumeric(rangeRelSD) && isscalar(rangeRelSD) && rangeRelSD >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for rangeRelSD! Needs to be a real non-negative scalar!');
            end
            if canBuildScenarioComponents(this.scenarioDimensionActive,this.numOfBeams)
                this.validateScenarioComponents(this.shiftSD,this.rangeAbsSD,rangeRelSD, ...
                    this.scenarioDimensionActive,this.numOfBeams,this.gantryAngleSD,this.couchAngleSD);
            end
            this.rangeRelSD = rangeRelSD;
            this.requestScenarioUpdate();
        end

        function set.rangeAbsSD(this,rangeAbsSD)
            valid = isnumeric(rangeAbsSD) && isscalar(rangeAbsSD) && rangeAbsSD >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for rangeAbsSD! Needs to be a real non-negative scalar!');
            end
            if canBuildScenarioComponents(this.scenarioDimensionActive,this.numOfBeams)
                this.validateScenarioComponents(this.shiftSD,rangeAbsSD,this.rangeRelSD, ...
                    this.scenarioDimensionActive,this.numOfBeams,this.gantryAngleSD,this.couchAngleSD);
            end
            this.rangeAbsSD = rangeAbsSD;
            this.requestScenarioUpdate();
        end

        function set.shiftSD(this,shiftSD)
            valid = isnumeric(shiftSD) && isrow(shiftSD) && numel(shiftSD) == 3 && all(shiftSD >= 0);
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for shiftSD! Needs to be 3-element non-negative numeric row vector!');
            end
            if canBuildScenarioComponents(this.scenarioDimensionActive,this.numOfBeams)
                this.validateScenarioComponents(shiftSD,this.rangeAbsSD,this.rangeRelSD, ...
                    this.scenarioDimensionActive,this.numOfBeams,this.gantryAngleSD,this.couchAngleSD);
            end
            this.shiftSD = shiftSD;
            this.requestScenarioUpdate();
        end

        function set.gantryAngleSD(this,gantryAngleSD)
            valid = isnumeric(gantryAngleSD) && isscalar(gantryAngleSD) && gantryAngleSD >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for gantryAngleSD! Needs to be a real non-negative scalar!');
            end
            if canBuildScenarioComponents(this.scenarioDimensionActive,this.numOfBeams)
                this.validateScenarioComponents(this.shiftSD,this.rangeAbsSD,this.rangeRelSD, ...
                    this.scenarioDimensionActive,this.numOfBeams,gantryAngleSD,this.couchAngleSD);
            end
            this.gantryAngleSD = gantryAngleSD;
            this.requestScenarioUpdate();
        end

        function set.couchAngleSD(this,couchAngleSD)
            valid = isnumeric(couchAngleSD) && isscalar(couchAngleSD) && couchAngleSD >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for couchAngleSD! Needs to be a real non-negative scalar!');
            end
            if canBuildScenarioComponents(this.scenarioDimensionActive,this.numOfBeams)
                this.validateScenarioComponents(this.shiftSD,this.rangeAbsSD,this.rangeRelSD, ...
                    this.scenarioDimensionActive,this.numOfBeams,this.gantryAngleSD,couchAngleSD);
            end
            this.couchAngleSD = couchAngleSD;
            this.requestScenarioUpdate();
        end

        function set.numOfBeams(this,numOfBeams)
            valid = isnumeric(numOfBeams) && isscalar(numOfBeams) && isfinite(numOfBeams) && ...
                round(numOfBeams) == numOfBeams && numOfBeams >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for numOfBeams! Needs to be a non-negative integer scalar!');
            end
            if canBuildScenarioComponents(this.scenarioDimensionActive,numOfBeams)
                this.validateScenarioComponents(this.shiftSD,this.rangeAbsSD,this.rangeRelSD, ...
                    this.scenarioDimensionActive,numOfBeams,this.gantryAngleSD,this.couchAngleSD);
            end
            this.numOfBeams = numOfBeams;
            if this.hasActiveAngularScenarioDimension()
                this.requestScenarioUpdate();
            end
        end

        function set.wcSigma(this,wcSigma)
            valid = isnumeric(wcSigma) && isscalar(wcSigma) && wcSigma >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for wcSigma! Needs to be a real positive scalar!');
            end
            this.wcSigma = wcSigma;
            this.requestScenarioUpdate();
        end

        function set.ctScenProb(this,ctScenProb)
            valid = isnumeric(ctScenProb) && ismatrix(ctScenProb) && ~isempty(ctScenProb) && ...
                size(ctScenProb,2) == 2 && ...
                all(round(ctScenProb(:,1)) == ctScenProb(:,1)) && all(ctScenProb(:,1) >= 0) && ...
                all(ctScenProb(:,2) > 0);
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError(['Invalid value for used ctScenProb! Needs to be a valid 2-column matrix ' ...
                    'with left column representing the scenario index and right column representing positive probabilities!']);
            end
            this.ctScenProb = ctScenProb;
            this.requestScenarioUpdate();
        end

        function set.scenarioDimensionActive(this,scenarioDimensionActive)
            scenarioDimensionActive = matRad_normalizeScenarioDimensionActive(scenarioDimensionActive);
            this.validateScenarioDimensionSupport(scenarioDimensionActive);
            if canBuildScenarioComponents(scenarioDimensionActive,this.numOfBeams)
                this.validateScenarioComponents(this.shiftSD,this.rangeAbsSD,this.rangeRelSD, ...
                    scenarioDimensionActive,this.numOfBeams,this.gantryAngleSD,this.couchAngleSD);
            end
            this.scenarioDimensionActive = scenarioDimensionActive;
            this.requestScenarioUpdate();
        end


        function scenarios = updateScenarios(this)
            %This function will always update the scenarios given the
            %current property settings

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('This abstract function needs to be implemented!');
        end

        function newInstance = extractSingleScenario(this,scenarioId)
            newInstance = matRad_NominalScenario();

            scenarioRowIx = this.resolveScenarioRowIx(scenarioId);
            ctScenId = this.getCtScenario(scenarioId);
            ctScenProbIx = find(this.ctScenProb(:,1) == ctScenId,1,'first');
            if isempty(ctScenProbIx)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Could not find CT scenario %d in ctScenProb.',ctScenId);
            end

            %First set properties that force an update
            newInstance.numOfCtScen         = 1;
            newInstance.numOfBeams          = this.numOfBeams;
            newInstance.ctScenProb          = this.ctScenProb(ctScenProbIx,:);
            newInstance.scenarioDimensionActive = this.scenarioDimensionActive;

            %Now overwrite existing variables for correct probabilties and
            %error realizations
            newInstance.numOfAvailableCtScen = this.numOfAvailableCtScen;
            singleMaskSize = [this.numOfAvailableCtScen ones(1,max(2,size(this.linearMask,2)-1))];
            scenMask = false(singleMaskSize);
            maskSubscripts = num2cell([ctScenId ones(1,numel(singleMaskSize)-1)]);
            scenMask(maskSubscripts{:}) = true;
            scenForProb = [ctScenId this.scenarioValues(scenarioRowIx,:)];
            newInstance.setScenarioRealizations(this.scenarioComponents,this.scenarioValues(scenarioRowIx,:), ...
                ctScenId,this.scenProb(scenarioRowIx),this.scenWeight(scenarioRowIx), ...
                scenForProb,[ctScenId ones(1,numel(singleMaskSize)-1)],scenMask);
            newInstance.maxAbsRangeShift = max(abs(this.getValue(scenarioId,'range.absolute')));
            newInstance.maxRelRangeShift = max(abs(this.getValue(scenarioId,'range.relative')));
            %newInstance.updateScenarios();
        end

        function scenIx = sub2scenIx(this,ctScenInput,setupScenIx,rangeScenIx,ctScenReference)
            warnDeprecatedScenarioApi('sub2scenIx');
            %Returns linear index in the scenario cell array from scenario
            %subscript indices. The optional ctScenReference disambiguates
            %whether ctScenInput is a local position ('position', default) or
            %an absolute CT scenario id ('id').
            if nargin < 5 || isempty(ctScenReference)
                ctScenReference = 'position';
            end
            ctScenId = resolveCtScenarioId(this,ctScenInput,ctScenReference);
            scenIx = scenarioSub2Ind(this,ctScenId,setupScenIx,rangeScenIx);
        end

        function scenarioRowIx = scenNum(this,fullScenIx)
            warnDeprecatedScenarioApi('scenNum');
            scenarioRowIx = this.getScenarioRowIndexFromDijIndex(fullScenIx);
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

        function scenRay = applyNonGeometricScenarioToRay(this,scenarioId,ray)
            ctScenId = this.getCtScenario(scenarioId);
            scenRay = ray;
            scenRay.radDepths = scenRay.radDepths{ctScenId};
            rangeApplicator = matRad_RangeShiftApplicator();
            scenRay.radDepths = rangeApplicator.applyToRadDepths(this,scenarioId,scenRay.radDepths);
            scenRay.radialDist_sq = scenRay.radialDist_sq{ctScenId};
            scenRay.ix = scenRay.ix{ctScenId};

            if isfield(scenRay,'geoDepths')
                scenRay.geoDepths = scenRay.geoDepths{ctScenId};
            end

            if isfield(scenRay,'latDists')
                scenRay.latDists = scenRay.latDists{ctScenId};
            end

            if isfield(scenRay,'isoLatDists')
                scenRay.isoLatDists = scenRay.isoLatDists{ctScenId};
            end
        end

        function ids = getNominalScenarioIds(this)
            ids = this.scenarioIdList(all(abs(this.scenarioValues) <= eps,2));
        end

        function fullScenIx = getDijScenarioIndex(this,scenarioId)
            scenarioRowIx = this.resolveScenarioRowIx(scenarioId);
            fullScenIx = find(this.scenMask);
            fullScenIx = fullScenIx(scenarioRowIx);
        end

        function fullScenIx = getDijScenarioIndexBySubscripts(this,ctScenInput,setupScenIx,rangeScenIx,varargin)
            [gantryScenIx,couchScenIx,ctScenReference,angularSubscriptsProvided] = ...
                parseScenarioSubscriptArgs(varargin{:});
            ctScenId = resolveCtScenarioId(this,ctScenInput,ctScenReference);
            fullScenIx = scenarioSub2Ind(this,ctScenId,setupScenIx,rangeScenIx, ...
                gantryScenIx,couchScenIx,angularSubscriptsProvided);
        end

        function tf = isScenarioActiveBySubscripts(this,ctScenInput,setupScenIx,rangeScenIx,varargin)
            fullScenIx = this.getDijScenarioIndexBySubscripts(ctScenInput,setupScenIx,rangeScenIx,varargin{:});
            tf = this.scenMask(fullScenIx);
        end

        function scenarioRowIx = getScenarioRowIndexFromDijIndex(this,fullScenIx)
            scenarioRowIx = find(find(this.scenMask) == fullScenIx,1,'first');
        end

        function scenarioRowIx = getScenarioNumberFromDijIndex(this,fullScenIx)
            warnDeprecatedScenarioApi('getScenarioNumberFromDijIndex');
            scenarioRowIx = this.getScenarioRowIndexFromDijIndex(fullScenIx);
        end

        function setupShift = getSetupShiftByIndex(this,setupScenIx)
            scenarioRowIx = find(this.linearMask(:,2) == setupScenIx,1,'first');
            if isempty(scenarioRowIx)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Setup scenario index %d does not exist.',setupScenIx);
            end
            setupShift = this.getSetupShift(this.scenarioIdList(scenarioRowIx));
        end

        function [absRangeShift,relRangeShift] = getRangeShiftByIndex(this,rangeScenIx)
            scenarioRowIx = find(this.linearMask(:,3) == rangeScenIx,1,'first');
            if isempty(scenarioRowIx)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Range scenario index %d does not exist.',rangeScenIx);
            end
            rangeShift = this.getRangeShift(this.scenarioIdList(scenarioRowIx));
            absRangeShift = rangeShift(1);
            relRangeShift = rangeShift(2);
        end

        function fp = fingerprint(this)
            data = struct();
            data.name = this.name;
            data.componentNames = this.scenarioValueNames;
            data.componentActive = [this.scenarioComponents.active];
            data.values = this.scenarioValues;
            data.ctScenIds = this.scenarioCtScenIds;
            data.weights = this.scenWeight;
            fp = stableScenarioFingerprint(data);
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

        function setScenarioRealizations(this,components_,scenarioValues_,ctScenIds_, ...
                scenProb_,scenWeight_,scenForProb_,linearMask_,scenMask_)
            this.scenarioComponents = components_;
            this.scenarioValueNames = {components_.name};
            this.scenarioValues = scenarioValues_;
            this.scenarioIdList = (1:size(scenarioValues_,1))';
            this.scenarioCtScenIds = ctScenIds_(:);
            this.scenarioApplicators = { ...
                matRad_CtScenarioApplicator(), ...
                matRad_SetupShiftApplicator(), ...
                matRad_RangeShiftApplicator(), ...
                matRad_GantryAngleApplicator(), ...
                matRad_CouchAngleApplicator()};

            this.isoShift = this.extractScenarioColumns({'setup.x','setup.y','setup.z'});
            this.absRangeShift = this.extractScenarioColumns({'range.absolute'});
            this.relRangeShift = this.extractScenarioColumns({'range.relative'});
            this.gantryAngleOffset = this.extractScenarioColumns(this.beamAngleComponentNames('gantry'));
            this.couchAngleOffset = this.extractScenarioColumns(this.beamAngleComponentNames('couch'));
            this.maxAbsRangeShift = max(this.absRangeShift);
            this.maxRelRangeShift = max(this.relRangeShift);
            this.scenForProb = scenForProb_;
            this.scenProb = scenProb_(:);
            this.scenWeight = scenWeight_(:);
            this.linearMask = linearMask_;
            this.scenMask = scenMask_;
        end

        %% Deprecated functions / properties
        function newInstance = extractSingleNomScen(this,~,scenIdx)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning(['The function extractSingleNomScen of the scenario class will soon be deprecated! ' ...
                'Use extractSingleScenario instead!']);
            newInstance = this.extractSingleScenario(scenIdx);
        end

        function t = TYPE(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property TYPE of the scenario class will soon be deprecated!');
            t = this.name;
        end

        function value = get.wcFactor(this)
            value = this.wcSigma;
        end

        function value = get.ctScenIx(this)
            warnDeprecatedScenarioApi('ctScenIx');
            value = this.scenarioCtScenIds;
        end

        function value = get.scenarioCtScen(this)
            warnDeprecatedScenarioApi('scenarioCtScen');
            value = this.scenarioCtScenIds;
        end

        function set.wcFactor(this,value)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property wcFactor of the scenario class will soon be deprecated!');
            this.wcSigma = value;
        end

    end

    methods (Access = protected)
        function initializeScenarioModel(this)
            this.updateScenariosOnChange = true;
            this.updateScenarios();
        end

        function requestScenarioUpdate(this)
            if this.updateScenariosOnChange && canBuildScenarioComponents(this.scenarioDimensionActive,this.numOfBeams)
                this.updateScenarios();
            end
        end

        function components = getScenarioComponents(this)
            components = matRad_createScenarioComponents(this.shiftSD,this.rangeAbsSD, ...
                this.rangeRelSD,this.scenarioDimensionActive,this.numOfBeams, ...
                this.gantryAngleSD,this.couchAngleSD);
        end

        function validateScenarioComponents(this,shiftSD,rangeAbsSD,rangeRelSD, ...
                scenarioDimensionActive,numOfBeams,gantryAngleSD,couchAngleSD)
            matRad_createScenarioComponents(shiftSD,rangeAbsSD,rangeRelSD, ...
                scenarioDimensionActive,numOfBeams,gantryAngleSD,couchAngleSD);
        end

        function validateScenarioDimensionSupport(~,~)
            % Subclasses can reject dimensions they do not implement.
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
        %{
        %TODO: implement automatic collection of available scenario classes

        function metaScenarioModels = getAvailableModels()
            matRad_cfg = MatRad_Config.instance();

            %Use the root folder and the scenarios folder only
            folders = {matRad_cfg.matRadRoot,mfilename("fullpath")};

            %
        end
        %}

        function types = AvailableScenCreationTYPE()
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The function/property AvailableScenarioCreationTYPE of the scenario class will soon be deprecated!');
            %Hardcoded for legacy scenario type compatibility
            types = {'nomScen','wcScen','impScen','truncatedImpScen','rndScen'};
        end
    end
end

function tf = canBuildScenarioComponents(scenarioDimensionActive,numOfBeams)

scenarioDimensionActive = matRad_normalizeScenarioDimensionActive(scenarioDimensionActive);
hasAngularDimension = any(strcmp(scenarioDimensionActive,'gantry')) || ...
    any(strcmp(scenarioDimensionActive,'couch'));
tf = ~hasAngularDimension || numOfBeams >= 1;

end

function ctScenId = resolveCtScenarioId(scenarioModel,ctScenInput,ctScenReference)

validatePositiveIntegerScalar(ctScenInput,'ctScen');
ctScenReference = normalizeCtScenReference(ctScenReference);

switch ctScenReference
    case 'position'
        ctScenIx = ctScenInput;
        if ctScenIx > size(scenarioModel.ctScenProb,1)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('CT scenario position %d exceeds the scenario model size.',ctScenIx);
        end
        ctScenId = scenarioModel.ctScenProb(ctScenIx,1);
    case 'id'
        ctScenId = ctScenInput;
        if ~any(scenarioModel.ctScenProb(:,1) == ctScenId)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Could not find CT scenario %d in the scenario model.',ctScenId);
        end
end

end

function ctScenReference = normalizeCtScenReference(ctScenReference)

if isstring(ctScenReference) && isscalar(ctScenReference)
    ctScenReference = char(ctScenReference);
end

if ~ischar(ctScenReference)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('ctScenReference must be ''position'' or ''id''.');
end

switch lower(ctScenReference)
    case {'position','ctscenposition'}
        ctScenReference = 'position';
    case {'id','ctscenid'}
        ctScenReference = 'id';
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('ctScenReference must be ''position'' or ''id''.');
end

end

function [gantryScenIx,couchScenIx,ctScenReference,angularSubscriptsProvided] = parseScenarioSubscriptArgs(varargin)

gantryScenIx = 1;
couchScenIx = 1;
ctScenReference = 'position';
angularSubscriptsProvided = false;

if nargin > 0
    if ischar(varargin{end}) || (isstring(varargin{end}) && isscalar(varargin{end}))
        ctScenReference = varargin{end};
        varargin(end) = [];
    end
end

if numel(varargin) >= 1 && ~isempty(varargin{1})
    gantryScenIx = varargin{1};
    angularSubscriptsProvided = true;
end

if numel(varargin) >= 2 && ~isempty(varargin{2})
    couchScenIx = varargin{2};
    angularSubscriptsProvided = true;
end

if numel(varargin) > 2
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Too many scenario subscript arguments.');
end

end

function scenIx = scenarioSub2Ind(scenarioModel,ctScenId,setupScenIx,rangeScenIx,gantryScenIx,couchScenIx,angularSubscriptsProvided)

if nargin < 5 || isempty(gantryScenIx)
    gantryScenIx = 1;
end

if nargin < 6 || isempty(couchScenIx)
    couchScenIx = 1;
end

if nargin < 7 || isempty(angularSubscriptsProvided)
    angularSubscriptsProvided = false;
end

validatePositiveIntegerScalar(setupScenIx,'setupScenIx');
validatePositiveIntegerScalar(rangeScenIx,'rangeScenIx');
validatePositiveIntegerScalar(gantryScenIx,'gantryScenIx');
validatePositiveIntegerScalar(couchScenIx,'couchScenIx');

scenMaskSize = size(scenarioModel.scenMask);
if numel(scenMaskSize) < 5
    scenMaskSize(end + 1:5) = 1;
end

if ctScenId > scenMaskSize(1) || setupScenIx > scenMaskSize(2) || ...
        rangeScenIx > scenMaskSize(3) || gantryScenIx > scenMaskSize(4) || ...
        couchScenIx > scenMaskSize(5)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Scenario subscript exceeds the scenario mask dimensions.');
end

if angularSubscriptsProvided || numel(size(scenarioModel.scenMask)) <= 3
    scenIx = sub2ind(scenMaskSize,ctScenId,setupScenIx,rangeScenIx,gantryScenIx,couchScenIx);
else
    activeSlice = scenarioModel.scenMask(ctScenId,setupScenIx,rangeScenIx,:,:);
    activeInSlice = find(activeSlice(:),1,'first');
    if isempty(activeInSlice)
        scenIx = sub2ind(scenMaskSize,ctScenId,setupScenIx,rangeScenIx,1,1);
    else
        [gantryScenIx,couchScenIx] = ind2sub([scenMaskSize(4) scenMaskSize(5)],activeInSlice);
        scenIx = sub2ind(scenMaskSize,ctScenId,setupScenIx,rangeScenIx,gantryScenIx,couchScenIx);
    end
end

end

function validatePositiveIntegerScalar(value,valueName)

if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
        round(value) == value && value >= 1)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('%s must be a positive integer scalar.',valueName);
end

end

function warnDeprecatedScenarioApi(symbol)

persistent warnedSymbols;

if isempty(warnedSymbols)
    warnedSymbols = {};
end

if any(strcmp(warnedSymbols,symbol))
    return;
end

warnedSymbols{end+1} = symbol;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispWarning(['The scenario function/property "%s" will be deprecated soon! ' ...
    '\nCheck out the new scenario realization API.'],symbol);

end

function fp = stableScenarioFingerprint(data)

try
    encoded = jsonencode(data);
catch
    encoded = evalc('disp(data)');
end

encoded = uint8(encoded);
hash = uint32(2166136261);
prime = uint32(16777619);
for i = 1:numel(encoded)
    hash = bitxor(hash,uint32(encoded(i)));
    hash = uint32(mod(double(hash) * double(prime),2^32));
end
fp = sprintf('%08x',hash);

end
