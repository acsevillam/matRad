function [cstStat, doseStat, meta, gammaFig, robustnessFig1, robustnessFig2] = matRad_samplingAnalysis(ct,cst,pln,caSampRes,mSampDose,resultGUInomScen,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad uncertainty sampling analysis function
% 
% call
%   [structureStat, doseStat] = samplingAnalysis(ct,cst,subIx,mSampDose,w)
%
% input
%   ct:                 ct cube
%   cst:                matRad cst struct
%   pln:                matRad's pln struct
%   caSampRes:          cell array of sampling results depicting plan 
%                       parameter
%   mSampDose:          matrix holding the sampled doses, each row 
%                       corresponds to one dose sample                      
%   resultGUInomScen:   resultGUI struct of the nominal plan
%   varargin:           optional Name/Value pairs for additional custom
%                       settings
%                       - 'phaseProb':      CT scenario phase probabilities
%                       - 'gammaCriterion': 1x2 vector [%  mm]
%                       - 'gammaCriteria':  legacy alias for gammaCriterion
%                       - 'robustnessCriteria': 1x2 vector [%  %]
%                       - 'slice':          CT slice used to create figures
%                       Dose statistics are evaluated per fraction.
%                       - 'Percentiles':    vector with desired percentiles
%                                           between (0,1)
%   
%
% output
%   cstStat         structure-wise statistics (mean, max, percentiles, ...)
%   doseStat        dose-wise statistics (mean, max, percentiles, ...)
%   meta            contains additional information about sampling analysis
%   gammaFig        gamma analysis figure if slice was provided
%   robustnessFig1  robustness analysis figure for combined deviation/std metric
%   robustnessFig2  robustness analysis figure for binary robustness criterion
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check integrity of statistics
matRad_cfg = MatRad_Config.instance();
gammaFig = [];
robustnessFig1 = [];
robustnessFig2 = [];

phaseProb = [];
if ~isempty(varargin) && isnumeric(varargin{1}) && isvector(varargin{1})
    phaseProb = varargin{1};
    varargin(1) = [];
end

p = inputParser;
p.CaseSensitive = false;
p.addParameter('gammaCriterion',[2 2],@(g) numel(g) == 2 && isnumeric(g) && all(g > 0));
p.addParameter('gammaCriteria',[2 2],@(g) numel(g) == 2 && isnumeric(g) && all(g > 0));
p.addParameter('robustnessCriteria',[5 5],@(r) numel(r) == 2 && isnumeric(r) && all(r > 0));
p.addParameter('phaseProb',phaseProb,@(phaseProbValue) isempty(phaseProbValue) || (isnumeric(phaseProbValue) && isvector(phaseProbValue) && all(phaseProbValue >= 0)));
p.addParameter('slice',[],@(slice) isempty(slice) || (isnumeric(slice) && isscalar(slice) && slice >= 0));
p.addParameter('percentiles',[0.01 0.05 0.125 0.25 0.5 0.75 0.875 0.95 0.99],@(p) (isscalar(p) || isvector(p)) && isnumeric(p) && all(p > 0 & p < 1)); 

parse(p,varargin{:});

meta = p.Results;
if any(strcmp(p.UsingDefaults,'gammaCriterion')) && ~any(strcmp(p.UsingDefaults,'gammaCriteria'))
    meta.gammaCriterion = meta.gammaCriteria;
end
meta.gammaCriteria = meta.gammaCriterion;
meta.doseMode = 'perFraction';

meta.sufficientStatistics = matRad_checkSampIntegrity(pln.multScen);
calcRobustness = nargout > 3 || ~any(strcmp(p.UsingDefaults,'robustnessCriteria'));

ix    = pln.subIx;
scenWeights = getScenarioWeights(pln,meta.phaseProb);

%% generate structurewise statistics
cstStat = struct();
% reassing dvh to stats structure
for i = 1:size(resultGUInomScen.cst,1)
    cstStat(i).name = caSampRes(1).qi(i).name;
    for l = 1:numel(caSampRes)
        % check if structures still match
        if any(~strcmp(cstStat(i).name, {caSampRes(l).dvh(i).name, caSampRes(l).qi(i).name}))
            matRad_cfg.dispError('matRad: Error, wrong structure.');
        end
        cstStat(i).dvh(l).doseGrid     = caSampRes(l).dvh(i).doseGrid;
        cstStat(i).dvh(l).volumePoints = caSampRes(l).dvh(i).volumePoints;
        cstStat(i).qi(l)               = caSampRes(l).qi(i);
        cstStat(i).w(l)                = scenWeights(l)';
    end  
end

if calcRobustness
    refDose = getTargetReferenceDose(cst,pln);
end

%% calculate mean and std cube
% compute doseMatrix with columns correspond to scenarios

doseStat.meanCube              = zeros(ct.cubeDim);
doseStat.stdCube               = zeros(ct.cubeDim);

doseStat.meanCubeW             = zeros(ct.cubeDim);
doseStat.stdCubeW              = zeros(ct.cubeDim);

doseStat.meanCube(ix)  = mean(mSampDose,2);
doseStat.stdCube(ix)   = std(mSampDose,1,2);
doseStat.meanCubeW(ix) = wMean(mSampDose',scenWeights)';
doseStat.stdCubeW(ix)  = wStd(mSampDose',scenWeights)';
    

% gamma cube
doseCube = resultGUInomScen.(pln.bioParam.quantityVis);
if strncmp(pln.bioParam.quantityVis,'RBExD', 5)
    doseStat.gammaAnalysis.cube1Name = 'resultGUInomScen.RBExD';
else
    doseStat.gammaAnalysis.cube1Name = 'resultGUInomScen.physicalDose';
end

doseStat.gammaAnalysis.cube1 = doseCube;
doseStat.gammaAnalysis.cube2 = doseStat.meanCubeW;
doseStat.gammaAnalysis.cube2Name = 'doseStat.meanCubeW';

matRad_cfg.dispInfo(['matRad: Performing gamma index analysis with parameters ', num2str(meta.gammaCriterion), ' [%% mm] \n']);
doseStat.gammaAnalysis.doseAgreement = meta.gammaCriterion(1);
doseStat.gammaAnalysis.distAgreement = meta.gammaCriterion(2);

[doseStat.gammaAnalysis.gammaCube,doseStat.gammaAnalysis.gammaPassRate,gammaFig] = matRad_gammaIndex(doseCube,doseStat.meanCubeW,[ct.resolution.x ct.resolution.y ct.resolution.z],meta.gammaCriterion,meta.slice);

if calcRobustness
    % robustness cube
    if strncmp(pln.bioParam.quantityVis,'RBExD', 5)
        doseStat.robustnessAnalysis.cubeName = 'resultGUInomScen.RBExD';
    else
        doseStat.robustnessAnalysis.cubeName = 'resultGUInomScen.physicalDose';
    end

    doseStat.robustnessAnalysis.meanCubeW = doseStat.meanCubeW;
    doseStat.robustnessAnalysis.stdCubeW = doseStat.stdCubeW;
    doseStat.robustnessAnalysis.refDose = refDose;
    doseStat.robustnessAnalysis.criteria = meta.robustnessCriteria;

    if isfinite(refDose)
        matRad_cfg.dispInfo(['matRad: Performing standard deviation analysis with parameters ', num2str(meta.robustnessCriteria), ' [%% %%] \n']);

        [doseStat.robustnessAnalysis.robustnessCube1,doseStat.robustnessAnalysis.robPassRate1,robustnessFig1] = calcRobustnessIndex1(doseStat.meanCubeW,doseStat.stdCubeW,refDose,meta.robustnessCriteria,meta.slice,ct,cst);
        doseStat.robustnessAnalysis.robustnessIndex1 = doseStat.robustnessAnalysis.robPassRate1/100;

        [doseStat.robustnessAnalysis.robustnessCube2,doseStat.robustnessAnalysis.robPassRate2,robustnessFig2] = calcRobustnessIndex2(doseStat.meanCubeW,doseStat.stdCubeW,refDose,meta.robustnessCriteria,meta.slice,ct,cst);
        doseStat.robustnessAnalysis.robustnessIndex2 = doseStat.robustnessAnalysis.robPassRate2/100;
    else
        matRad_cfg.dispWarning('Could not determine target reference dose. Skipping robustness index calculation.\n');
        doseStat.robustnessAnalysis.robustnessCube1 = [];
        doseStat.robustnessAnalysis.robPassRate1 = [];
        doseStat.robustnessAnalysis.robustnessIndex1 = [];
        doseStat.robustnessAnalysis.robustnessCube2 = [];
        doseStat.robustnessAnalysis.robPassRate2 = [];
        doseStat.robustnessAnalysis.robustnessIndex2 = [];
    end
end

%% percentiles
percentiles     = meta.percentiles;
percentileNames = cell(numel(percentiles),1);
% create fieldnames
for i = 1:numel(percentiles)
    percentileNames{i} = ['P',num2str(percentiles(i)*100)];
end
% create table rownames
metric = vertcat({'mean';'min';'max';'std'},percentileNames{:});

% create statstics where structure based results (QI and DVH) are available
for i = 1:size(cst,1)
    cstStat(i).percentiles = percentiles;
    cstStat(i).metric      = metric;
    cstStat(i).dvhStat     = calcDVHStat(cstStat(i).dvh,cstStat(i).percentiles,cstStat(i).w);
    cstStat(i).qiStat      = calcQiStat(cstStat(i).qi,cstStat(i).percentiles,cstStat(i).w);
end


% dvh statistics
    function [dvhStat, doseGrid] = calcDVHStat(dvh,percentiles,w)
        doseGrid = dvh(1).doseGrid;
        dvhMat = NaN * ones(numel(dvh),numel(dvh(1).volumePoints));
        for j = 1:numel(dvh)
            dvhMat(j,:) = dvh(j).volumePoints;            
        end
        % for statistical reasons, treat NaN as 0
        dvhMat(isnan(dvhMat)) = 0;
        
        dvhStat.mean.doseGrid     = doseGrid;
        dvhStat.mean.volumePoints = wMean(dvhMat,w);
        % [~,argmin] = min(dvhStat.mean.volumePoints);
        % dvhStat.mean(2,argmin + 1:end) = NaN;
        
        dvhStat.min.doseGrid     = doseGrid;
        dvhStat.min.volumePoints = min(dvhMat);
        % [~,argmin] = min(dvhStat.min.volumePoints);
        % dvhStat.min(2,argmin + 1:end) = NaN;
        
        dvhStat.max.doseGrid     = doseGrid;
        dvhStat.max.volumePoints = max(dvhMat);
        % [~,argmin] = min(dvhStat.max.volumePoints);
        % dvhStat.max(2,argmin + 1:end) = NaN;
        
        dvhStat.std.doseGrid     = doseGrid;
        dvhStat.std.volumePoints = wStd(dvhMat,w);

        dvhStat.percDVH = NaN * ones(numel(percentiles),numel(doseGrid));
        
        for j = 1:size(dvhMat,2)
            wQ =  matRad_weightedQuantile(dvhMat(:,j), percentiles, w', false);
            dvhStat.percDVH(:,j) = wQ;
        end

    end % eof calcDVHStat

    % qi statistics
    function qiStat = calcQiStat(qi,percentiles,w)
        fields = fieldnames(qi);
        % remove name field
        if sum(strcmp('VOIname', fields)) >= 1
            qi = rmfield(qi, 'VOIname');
        end
        fields = fieldnames(qi);
        qiStruct = qi;
        
        % create helper matlab structure which will be converted to table
        qiStatH = struct();
        for j = 1:numel(fields)
            if numel([qiStruct(:).(fields{j})]) == numel(w)
                qiStatH(1).(fields{j}) = wMean([qiStruct(:).(fields{j})],w);
                qiStatH(2).(fields{j}) = min([qiStruct(:).(fields{j})]);
                qiStatH(3).(fields{j}) = max([qiStruct(:).(fields{j})]);
                qiStatH(4).(fields{j}) = wStd([qiStruct(:).(fields{j})],w);
                wQ = matRad_weightedQuantile([qiStruct(:).(fields{j})], percentiles, w', false);
                for k = 1:numel(wQ)
                    sIx = k + 4;
                    qiStatH(sIx).(fields{j}) = wQ(k);
                end
            else
                for k = 1:(4 + numel(percentiles))
                    qiStatH(k).(fields{j}) = [];
                end
            end
        end
        env = matRad_getEnvironment();
        if strcmp(env,'MATLAB')
            qiStat = struct2table(qiStatH);
            qiStat.Properties.RowNames = metric;
        else
            qiStat = qiStatH;
        end
    end % eof calcQiStat

    function S = wMean(X,w)
        if exist('w','var') && ~isempty(w)
            if isvector(X) && isvector(w)
                S = reshape(w,1,[]) * reshape(X,[],1) / sum(w);
            else
                % row-wise
                S = reshape(w,1,[]) * X ./ sum(w);        
            end

        else
            S = mean(X);
        end
    end

    function S = wStd(X,w)
        if exist('w','var') && ~isempty(w)
            if isvector(X) && isvector(w)
                x = reshape(X,[],1);
                w = reshape(w,[],1);
                mu = sum(w.*x)./sum(w);
                S = sqrt(sum(w.*(x - mu).^2)./sum(w));
            else
                mu = wMean(X,w);
                S = sqrt(reshape(w,1,[]) * (bsxfun(@minus,X,mu).^2) ./ sum(w));
            end
        else
            S = std(X,1);
        end
    end

    function scenWeights = getScenarioWeights(pln,phaseProb)
        hasScenWeight = (isobject(pln.multScen) && isprop(pln.multScen,'scenWeight')) || ...
            (isstruct(pln.multScen) && isfield(pln.multScen,'scenWeight'));

        if hasScenWeight && ~isempty(pln.multScen.scenWeight)
            scenWeights = pln.multScen.scenWeight(:);
        else
            scenWeights = pln.multScen.scenProb(:);
        end

        if numel(scenWeights) ~= numel(caSampRes)
            matRad_cfg.dispError('Number of scenario weights does not match number of sampled scenarios.');
        end

        if isempty(phaseProb)
            scenWeights = normalizeScenarioWeights(scenWeights);
            return;
        end

        phaseProb = phaseProb(:);

        if isempty(pln.multScen.linearMask)
            matRad_cfg.dispWarning('Phase probabilities were provided, but multScen.linearMask is empty. Using multScen.scenWeight.\n');
            scenWeights = normalizeScenarioWeights(scenWeights);
            return;
        end

        ctScen = pln.multScen.linearMask(:,1);
        if numel(phaseProb) < max(ctScen)
            matRad_cfg.dispError('Phase probability vector must contain one entry for every CT scenario.');
        end

        currentPhaseProb = ones(size(scenWeights));
        hasCtScenProb = (isobject(pln.multScen) && isprop(pln.multScen,'ctScenProb')) || ...
            (isstruct(pln.multScen) && isfield(pln.multScen,'ctScenProb'));
        if hasCtScenProb && ~isempty(pln.multScen.ctScenProb)
            for s = 1:numel(ctScen)
                phaseIx = find(pln.multScen.ctScenProb(:,1) == ctScen(s),1,'first');
                if ~isempty(phaseIx)
                    currentPhaseProb(s) = pln.multScen.ctScenProb(phaseIx,2);
                end
            end
        end

        invalidPhaseProb = currentPhaseProb <= 0;
        if any(invalidPhaseProb)
            matRad_cfg.dispWarning('Current CT scenario probability contains zeros. Using existing scenario weights for affected scenarios.\n');
            currentPhaseProb(invalidPhaseProb) = 1;
        end

        scenWeights = scenWeights ./ currentPhaseProb .* phaseProb(ctScen);
        scenWeights = normalizeScenarioWeights(scenWeights);
    end

    function scenWeights = normalizeScenarioWeights(scenWeights)
        scenWeights = scenWeights(:);
        scenWeights(~isfinite(scenWeights) | scenWeights < 0) = 0;
        if sum(scenWeights) <= 0
            matRad_cfg.dispError('Scenario weights must contain at least one positive finite value.');
        end
        scenWeights = scenWeights./sum(scenWeights);
    end

    function refDose = getTargetReferenceDose(cst,pln)
        refDose = inf;

        for i = 1:size(cst,1)
            if ~isequal(cst{i,3},'TARGET')
                continue;
            end

            objectives = cst{i,6};
            if isempty(objectives)
                continue;
            end

            if isstruct(objectives)
                objectives = num2cell(arrayfun(@matRad_DoseOptimizationFunction.convertOldOptimizationStruct,objectives));
            end

            for runObjective = 1:numel(objectives)
                obj = objectives{runObjective};
                if ~isa(obj,'matRad_DoseOptimizationFunction')
                    try
                        obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
                    catch ME
                        matRad_cfg.dispWarning('Objective/Constraint not valid!\n%s',ME.message);
                        continue;
                    end
                end

                isUnderdoseObjective = isa(obj,'DoseObjectives.matRad_SquaredDeviation') || ...
                    isa(obj,'DoseObjectives.matRad_SquaredUnderdosing') || ...
                    isa(obj,'DoseObjectives.matRad_SquaredBertoluzzaDeviation2') || ...
                    isa(obj,'DoseObjectives.matRad_MinDVH');

                if isUnderdoseObjective
                    doseParameters = obj.getDoseParameters();
                    doseParameters = doseParameters(isfinite(doseParameters));
                    if ~isempty(doseParameters)
                        refDose = min([refDose; doseParameters(:)]);
                    end
                end
            end
        end

        if isfinite(refDose)
            refDose = refDose/pln.numOfFractions;
        else
            matRad_cfg.dispWarning('Target has no objective that penalizes underdosage. Reference dose unavailable.\n');
        end
    end

    function [robCube,robPassRate,robustnessFig] = calcRobustnessIndex1(meanCube,stdCube,refDose,criteria,slice,ct,cst)
        meanDoseThreshold = criteria(1);
        stdThreshold = criteria(2);

        [targetMask,doseMask,refScen] = getRobustnessMasks(meanCube,stdCube,ct,cst);

        meanDoseCrit = abs(meanCube - refDose)/refDose*100/meanDoseThreshold;
        stdCrit = stdCube/refDose*100/stdThreshold;
        robCube = sqrt(meanDoseCrit.^2 + stdCrit.^2).*doseMask;

        targetIx = ~isnan(targetMask);
        numOfTargetVoxels = sum(targetIx(:));
        if numOfTargetVoxels > 0
            numOfPassRobustness = sum((robCube(:) <= 1) & targetIx(:));
            robPassRate = 100 * numOfPassRobustness / numOfTargetVoxels;
        else
            robPassRate = NaN;
        end

        robustnessFig = plotRobustnessIndex1(robCube,targetMask,doseMask,robPassRate,meanDoseThreshold,stdThreshold,slice,ct,cst,refScen);
    end

    function [robCube,robPassRate,robustnessFig] = calcRobustnessIndex2(meanCube,stdCube,refDose,criteria,slice,ct,cst)
        meanDoseThreshold = criteria(1);
        stdThreshold = criteria(2);

        [targetMask,doseMask,refScen] = getRobustnessMasks(meanCube,stdCube,ct,cst);

        meanDoseCrit = meanCube >= refDose;
        stdCrit = stdCube <= stdThreshold;
        robCube = meanDoseCrit & stdCrit;

        targetIx = ~isnan(targetMask);
        numOfTargetVoxels = sum(targetIx(:));
        if numOfTargetVoxels > 0
            numOfPassRobustness = sum((robCube(:) == 1) & targetIx(:));
            robPassRate = 100 * numOfPassRobustness / numOfTargetVoxels;
        else
            robPassRate = NaN;
        end

        robustnessFig = plotRobustnessIndex2(robCube,doseMask,robPassRate,meanDoseThreshold,stdThreshold,slice,ct,cst,refScen);
    end

    function [targetMask,doseMask,refScen] = getRobustnessMasks(meanCube,stdCube,ct,cst)
        if ~isfield(ct,'refScen')
            refScen = 1;
        else
            refScen = ct.refScen;
        end

        targetMask = NaN(size(meanCube));
        for i = 1:size(cst,1)
            if isequal(cst{i,3},'TARGET') && numel(cst{i,4}) >= refScen && ~isempty(cst{i,4}{refScen})
                targetMask(cst{i,4}{refScen}) = 1;
            end
        end

        doseMask = meanCube > 0 | stdCube > 0;
    end

    function robustnessFig = plotRobustnessIndex1(robCube,targetMask,doseMask,robPassRate,meanDoseThreshold,stdThreshold,slice,ct,cst,refScen)
        robustnessFig = [];
        if isempty(slice)
            return;
        end

        plane = 3;
        maxRob = 5.01;
        doseWindow = [0 maxRob];

        mMap1 = round(1/maxRob*256);
        mMap2 = 256 - mMap1;

        colormap1 = [linspace(0.40,1,mMap1)',linspace(1,1,mMap1)',linspace(0.40,1,mMap1)'];
        colormap2 = matRad_getColormap('gammaIndex',2*mMap2);
        myColormap = [colormap1; colormap2(mMap2+1:end-1,:)];

        robustnessFig = figure;
        robustnessFig.Position(3:4) = [800 400];

        subplot(1,2,1);
        numSlices = ct.cubeDim(3);
        matRad_plotSliceWrapper(gca,ct,cst,refScen,robCube.*targetMask,plane,slice,[],[],colorcube,myColormap,doseWindow,[],[],'Deviation-Uncertainty',[],'LineWidth',1.5);
        title({[num2str(robPassRate,5) '% of points pass robustness criterion (' ...
            num2str(meanDoseThreshold) '% / ' num2str(stdThreshold) '%)']});

        ax1 = gca;
        sliderStep = getSliderStep(numSlices);
        b = uicontrol('Parent',robustnessFig,'Style','slider','Position',[70,5,280,23], ...
            'value',slice,'min',1,'max',numSlices,'SliderStep',sliderStep);
        b.Callback = @(es,ed) matRad_plotSliceWrapper(ax1,ct,cst,refScen,robCube.*targetMask,plane,round(es.Value),[],[],colorcube,myColormap,doseWindow,[],[],'Deviation-Uncertainty',[],'LineWidth',1.5);

        subplot(1,2,2);
        matRad_plotSliceWrapper(gca,ct,cst,refScen,(robCube <= 1).*doseMask,plane,slice,[],[],colorcube,[],[0 2.01],[],[],[],[],'LineWidth',1.5);
        title({[num2str(robPassRate,5) '% of points pass robustness criterion (' ...
            num2str(meanDoseThreshold) '% / ' num2str(stdThreshold) '%)']});

        ax2 = gca;
        b = uicontrol('Parent',robustnessFig,'Style','slider','Position',[430,5,280,23], ...
            'value',slice,'min',1,'max',numSlices,'SliderStep',sliderStep);
        b.Callback = @(es,ed) matRad_plotSliceWrapper(ax2,ct,cst,refScen,(robCube <= 1).*doseMask,plane,round(es.Value),[],[],colorcube,[],[0 2.01],[],[],[],[],'LineWidth',1.5);
    end

    function robustnessFig = plotRobustnessIndex2(robCube,doseMask,robPassRate,meanDoseThreshold,stdThreshold,slice,ct,cst,refScen)
        robustnessFig = [];
        if isempty(slice)
            return;
        end

        plane = 3;
        maxRob = 2.01;
        doseWindow = [0 maxRob];

        robustnessFig = figure;
        robustnessFig.Position(3:4) = [400 400];

        numSlices = ct.cubeDim(3);
        matRad_plotSliceWrapper(gca,ct,cst,refScen,robCube.*doseMask,plane,slice,[],[],colorcube,[],doseWindow,[],[],'Deviation-Uncertainty',[],'LineWidth',1.5);
        title({[num2str(robPassRate,5) '% of points pass robustness criterion (' ...
            num2str(meanDoseThreshold) '% / ' num2str(stdThreshold) '%)']});

        ax1 = gca;
        sliderStep = getSliderStep(numSlices);
        b = uicontrol('Parent',robustnessFig,'Style','slider','Position',[70,5,280,23], ...
            'value',slice,'min',1,'max',numSlices,'SliderStep',sliderStep);
        b.Callback = @(es,ed) matRad_plotSliceWrapper(ax1,ct,cst,refScen,robCube.*doseMask,plane,round(es.Value),[],[],colorcube,[],doseWindow,[],[],'Deviation-Uncertainty',[],'LineWidth',1.5);
    end

    function sliderStep = getSliderStep(numSlices)
        if numSlices > 1
            sliderStep = [1/(numSlices-1) 1/(numSlices-1)];
        else
            sliderStep = [1 1];
        end
    end

    % check integrity of scenario analysis (i.e. check number of scenarios)
    function statCheck = matRad_checkSampIntegrity(multScen)
       
        if multScen.totNumScen > 20
            totalNum = true;
        else
            totalNum = false;
        end
        
        statCheck = totalNum; % * .... *
        
    end

end
