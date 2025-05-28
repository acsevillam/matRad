function [caSampRes, mSampDose, pln, resultGUInomScen, resultGUIsampledScen]  = matRad_sampling(ct,stf,cst,pln,w,structSel,multScen)
% matRad_randomSampling enables sampling multiple treatment scenarios
%
% call
%   [caSampRes,mSampDose,pln,resultGUInomScen,resultGUIsampledScen] = matRad_Sampling2(ct,cst,pln)
%
% input
%   ct:         matRad ct information struct
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   cst:        matRad cst struct
%   w:          optional (if no weights available in stf): bixel weight
%               vector
% output
%   caSampRes:              cell array of sampling results depicting plan parameter
%   mSampDose:              matrix holding the sampled doses, each row corresponds to
%                           one dose sample
%   pln:                    matRad pln struct containing sampling information
%   resultGUInomScen:       resultGUI struct of the nominal scenario
%   resultGUIsampledScen:   resultGUI struct for all sampled scenarios
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

% save nonSampling pln for nominal scenario calculation and add dummy fields
plnNominal = pln;
% create nominal scenario
plnNominal.multScen = matRad_multScen(ct,'nomScen');

ctSamp = ct;

% either use existing multScen struct or create new one
if exist('multScen','var') && ~isempty(multScen)
    pln.multScen = multScen;
else
    % create random scenarios for sampling
    pln.multScen = matRad_multScen(ctSamp,'rndScen'); % 'impSamp' or 'wcSamp'
    pln.multScen.numOfShiftScen = matRad_cfg.defaults.samplingScenarios * ones(3,1);
    pln.multScen.numOfRangeShiftScen = matRad_cfg.defaults.samplingScenarios;
end

matRad_cfg.dispInfo('Using %d samples in total \n',pln.multScen.totNumScen);

if ~isfield(ct,'refScen')
    refScen=1;
else
    refScen=ct.refScen;
end

% default ct scenario for sampling
matRad_cfg.dispInfo('Sampling will be performed on ct scenario: %d \n',refScen);

for  i = 1:size(cst,1)
    if isequal(cst{i,3},'TARGET')
        
        % loop over target objectives and get the lowest dose objective
        refDose = inf;
        
        if isstruct(cst{i,6})
            cst{i,6} = num2cell(arrayfun(@matRad_DoseOptimizationFunction.convertOldOptimizationStruct,cst{i,6}));
        end
        
        for runObjective = 1:numel(cst{i,6})
            % check if this is an objective that penalizes underdosing
            obj = cst{i,6}{runObjective};
            if ~isa(obj,'matRad_DoseOptimizationFunction')
                try
                    obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
                catch ME
                    matRad_cfg.dispWarning('Objective/Constraint not valid!\n%s',ME.message)
                    continue;
                end
            end
            
            if isa(obj,'DoseObjectives.matRad_SquaredDeviation') || isa(obj,'DoseObjectives.matRad_SquaredUnderdosing') || isa(obj,'DoseObjectives.matRad_MinDVH') || isa(obj,'DoseObjectives.matRad_SquaredBertoluzzaDeviation2')
                refDose = (min(obj.getDoseParameters(),refDose));%/pln.numOfFractions;
            end
        end
        
        if refDose == inf
            sprintf('%s%s',i,'Warning: target has no objective that penalizes underdosage, ');
        end
        
    end
end

V = [];
% define voxels for sampling
if ~exist('structSel', 'var') || sum(size(structSel)) == 0
    V = [cst{:,4}];
    % final voxel subset for sampling
    subIx = unique(vertcat(V{:}));
else
    for i=1:size(cst,1)
        for j = 1:numel(structSel)
            if strcmp(structSel{j}, cst{i,2})
                V = [V cst{i,4}{1}];
            end
        end
    end
    % final voxel subset for sampling
    subIx = V;
end

% disable structures for DVH plotting which are not completely in subIx
for i = 1:size(cst,1)
    if ~all(ismember(cst{i,4}{1},subIx))
        cst{i,5}.Visible = false;
    end
end

% define variable for storing scenario doses
mSampDose   = single(zeros(numel(subIx),pln.multScen.totNumScen,1));
StorageInfo = whos('mSampDose');
matRad_cfg.dispInfo('matRad: Realizations variable will need: %f GB \n',StorageInfo.bytes/1e9);

% check if parallel toolbox is installed and license can be checked out
try
    [FlagParallToolBoxLicensed,msg]  = license('checkout','Distrib_Computing_Toolbox');
    if ~FlagParallToolBoxLicensed
        matRad_cfg.dispWarning('Could not check out parallel computing toolbox. \n');
    end
    
catch
    FlagParallToolBoxLicensed  = false;
end

%% calculate nominal scenario
nomScenTimer     = tic;
resultGUInomScen = matRad_calcDoseDirect(ct,stf,plnNominal,cst,w);
nomScenTime      = toc(nomScenTimer);
matRad_cfg.dispInfo('Finished nominal Scenario Calculation. Computation time: %f h \n',round(nomScenTime / 3600));

refVol = [2 5 50 95 98];
refGy = floor(linspace(0,refDose,6)*10)/10;

resultGUInomScen.dvh = matRad_calcDVH(cst,resultGUInomScen.(pln.bioParam.quantityVis),'cum');
dvhPoints            = resultGUInomScen.dvh(1).doseGrid;
nomQi                = matRad_calcQualityIndicators(cst,pln,resultGUInomScen.(pln.bioParam.quantityVis)*pln.numOfFractions,refGy,refVol);

resultGUInomScen.qi  = nomQi;
resultGUInomScen.cst = cst;

%% perform parallel sampling
if FlagParallToolBoxLicensed

    nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
    
    % Fallback para pruebas locales
    if isnan(nWorkers) || nWorkers <= 0
        nCores = feature('numcores');
        nWorkers = max(1, nCores - 2);
    end

    fprintf('[matRad_sampling] SLURM_CPUS_PER_TASK: %d\n', nWorkers);

    p = gcp('nocreate');
    if ~isempty(p)
        if p.NumWorkers ~= nWorkers
            fprintf('[matRad_sampling] closing existing parpool with %d workers\n', p.NumWorkers);
            delete(p);
        end
    end

    % Intentar abrir el nuevo parpool
    if isempty(gcp('nocreate'))
        try
            parpool('local', nWorkers);
        catch ME
            warning('[matRad_sampling] Error at initialize parpool: %s', ME.message);
            rethrow(ME);
        end
    end

    % rough estimate of total computation time
    totCompTime = ceil(pln.multScen.totNumScen / nWorkers) * nomScenTime * 1.35;
    matRad_cfg.dispInfo(['Approximate Total calculation time: ', num2str(round(totCompTime / 3600)), ...
        'h. Estimated finish: ', datestr(datetime('now') + seconds(totCompTime)), '\n']);
    
    if exist('parfor_progress', 'file') == 2
        FlagParforProgressDisp = true;
        parfor_progress(pln.multScen.totNumScen);  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
    else
        matRad_cfg.dispInfo('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n');
        FlagParforProgressDisp = false;
    end
    
    resultGUIsampledScen.(pln.bioParam.quantityVis) = cell(size(pln.multScen.scenMask));
    
    for i = 1:pln.multScen.totNumScen
        
        [ctScen,shiftScen,RangeScen] = deal(pln.multScen.linearMask(i,1),pln.multScen.linearMask(i,2),pln.multScen.linearMask(i,3));
        
        shiftScenMask = find(squeeze(pln.multScen.scenMask(1,:,:)));
        indProb = sub2ind([pln.multScen.totNumShiftScen pln.multScen.totNumRangeScen],shiftScen,RangeScen);
        
        matRad_cfg.dispInfo([' \n Sampling scenario ', num2str(i), ...
            ' of ', num2str(pln.multScen.totNumScen), ' \n  \n']);
        
        % create nominal scenario
        plnSamp             = pln;
        plnSamp.multScen	= pln.multScen.extractSingleNomScen(1,find(shiftScenMask==indProb));
        ctSamp              = ct;
        ctSamp.numOfCtScen  = 1;
        ctSamp.cubeHU       = [];
        ctSamp.cubeHU{1}    = ct.cubeHU{ctScen};
        if isfield(ct,'dvf')
            ctSamp.dvf          = [];
            ctSamp.dvf{1}       = ct.dvf{ctScen};
        end
        cstSamp             = cst;

       [numOfStruct, ~] = size(cstSamp);
       for structure = 1:numOfStruct
            cstSamp{structure,4}=[];
            cstSamp{structure,4}=cst{structure,4}(refScen);
       end

        resultSamp                 = matRad_calcDoseDirect(ctSamp,stf,plnSamp,cstSamp,w);
        
        if isfield(ctSamp,'dvf')

            if ~isequal(size(resultSamp.(plnSamp.bioParam.quantityVis)),ctSamp.dvfDim) || ~isequal(ctSamp.dvfType,'pull')
                matRad_cfg.dispWarning('Dose cube and deformation vector field dimensions are not equal. \n');

                % Instantiate elastic registration
                metadata.nItera = 100;
                metadata.dvfType = 'pull';
                register = matRad_ElasticImageRegistration(ct,cst,1,metadata);
                clear 'metadata';
                
                % Calculate deformation vector field
                [ct] = register.calcDVF();
                ctSamp.dvfDim       = ct.dvfDim;
                ctSamp.dvf{1}       = ct.dvf{ctScen};
            end

            resultSamp.([plnSamp.bioParam.quantityVis '_deformed']) = imwarp(resultSamp.(plnSamp.bioParam.quantityVis), permute(ctSamp.dvf{1},[2 3 4 1]));

        end

        resultGUIsampledScen.(plnSamp.bioParam.quantityVis){ctScen,shiftScen,RangeScen} = resultSamp.(plnSamp.bioParam.quantityVis) ;
        

        if isfield(resultSamp,[plnSamp.bioParam.quantityVis '_deformed'])
            resultGUIsampledScen.([plnSamp.bioParam.quantityVis '_deformed']){ctScen,shiftScen,RangeScen} = resultSamp.([plnSamp.bioParam.quantityVis '_deformed']);
            sampledDose                = resultSamp.([plnSamp.bioParam.quantityVis '_deformed'])(subIx);
        else
            sampledDose                = resultSamp.(plnSamp.bioParam.quantityVis)(subIx);
        end

        mSampDose(:,i)             = single(reshape(sampledDose,[],1));
        caSampRes(i).bioParam      = plnSamp.bioParam;
        caSampRes(i).relRangeShift = plnSamp.multScen.relRangeShift;
        caSampRes(i).absRangeShift = plnSamp.multScen.absRangeShift;
        caSampRes(i).isoShift      = plnSamp.multScen.isoShift;

        if isfield(resultSamp,[plnSamp.bioParam.quantityVis '_deformed'])
            caSampRes(i).dvh = matRad_calcDVH(cstSamp,resultSamp.([plnSamp.bioParam.quantityVis '_deformed']),'cum',dvhPoints);
            caSampRes(i).qi  = matRad_calcQualityIndicators(cstSamp,plnSamp,resultSamp.([plnSamp.bioParam.quantityVis '_deformed'])*plnSamp.numOfFractions,refGy,refVol);
        else
            caSampRes(i).dvh = matRad_calcDVH(cstSamp,resultSamp.(plnSamp.bioParam.quantityVis),'cum',dvhPoints);
            caSampRes(i).qi  = matRad_calcQualityIndicators(cstSamp,plnSamp,resultSamp.(plnSamp.bioParam.quantityVis)*plnSamp.numOfFractions,refGy,refVol);
        end
        
    end
   
    if FlagParforProgressDisp
        parfor_progress(0);
    end

    if ~isempty(p)
        if p.NumWorkers ~= nWorkers
            fprintf('[matRad_sampling] closing existing parpool with %d workers\n', p.NumWorkers);
            delete(p);
        end
    end
    
else
    %% perform seriel sampling
    % rough estimate of total computation time
    totCompTime = size(pln.multScen.scenForProb,1) * nomScenTime * 1.1;
    try
        matRad_cfg.dispInfo(['Approximate Total calculation time: ', num2str(round(totCompTime / 3600)), ...
            'h. Estimated finish: ', datestr(datetime('now') + seconds(totCompTime)), '\n']);
    catch
        matRad_cfg.dispInfo(['Approximate Total calculation time: ', num2str(round(totCompTime / 3600)), '\n']);
    end
    
    resultGUIsampledScen.(pln.bioParam.quantityVis) = cell(size(pln.multScen.scenMask));
    
    for i = 1:pln.multScen.totNumScen
        
        [ctScen,shiftScen,RangeScen] = deal(pln.multScen.linearMask(i,1),pln.multScen.linearMask(i,2),pln.multScen.linearMask(i,3));
        
        shiftScenMask = find(squeeze(pln.multScen.scenMask(1,:,:)));
        indProb = sub2ind([pln.multScen.totNumShiftScen pln.multScen.totNumRangeScen],shiftScen,RangeScen);
        
        matRad_cfg.dispInfo([' \n Sampling scenario ', num2str(i), ...
            ' of ', num2str(pln.multScen.totNumScen), ' \n  \n']);
        
        % create nominal scenario
        plnSamp             = pln;
        plnSamp.multScen	= pln.multScen.extractSingleNomScen(1,find(shiftScenMask==indProb));
        ctSamp              = ct;
        ctSamp.numOfCtScen  = 1;
        ctSamp.cubeHU       = [];
        ctSamp.cubeHU{1}    = ct.cubeHU{ctScen};
        ctSamp.dvf          = [];
        ctSamp.dvf{1}       = ct.dvf{ctScen};
        cstSamp             = cst;

       [numOfStruct, ~] = size(cstSamp);
       for structure = 1:numOfStruct
            cstSamp{structure,4}=[];
            cstSamp{structure,4}=cst{structure,4}(refScen);
       end

        resultSamp                 = matRad_calcDoseDirect(ctSamp,stf,plnSamp,cstSamp,w);
        
        if isfield(ctSamp,'dvf')

            if ~isequal(size(resultSamp.(plnSamp.bioParam.quantityVis)),ctSamp.dvfDim) || ~isequal(ctSamp.dvfType,'pull')
                matRad_cfg.dispWarning('Dose cube and deformation vector field dimensions are not equal. \n');

                % Instantiate elastic registration
                metadata.nItera = 100;
                metadata.dvfType = 'pull';
                register = matRad_ElasticImageRegistration(ct,cst,1,metadata);
                clear 'metadata';
                
                % Calculate deformation vector field
                [ct] = register.calcDVF();
                ctSamp.dvfDim       = ct.dvfDim;
                ctSamp.dvf{1}       = ct.dvf{ctScen};
            end

            resultSamp.([plnSamp.bioParam.quantityVis '_deformed']) = imwarp(resultSamp.(plnSamp.bioParam.quantityVis), permute(ctSamp.dvf{1},[2 3 4 1]));

        end

        resultGUIsampledScen.(plnSamp.bioParam.quantityVis){ctScen,shiftScen,RangeScen} = resultSamp.(plnSamp.bioParam.quantityVis) ;
        

        if isfield(resultSamp,[plnSamp.bioParam.quantityVis '_deformed'])
            resultGUIsampledScen.([plnSamp.bioParam.quantityVis '_deformed']){ctScen,shiftScen,RangeScen} = resultSamp.([plnSamp.bioParam.quantityVis '_deformed']);
            sampledDose                = resultSamp.([plnSamp.bioParam.quantityVis '_deformed'])(subIx);
        else
            sampledDose                = resultSamp.(plnSamp.bioParam.quantityVis)(subIx);
        end

        mSampDose(:,i)             = single(reshape(sampledDose,[],1));
        caSampRes(i).bioParam      = plnSamp.bioParam;
        caSampRes(i).relRangeShift = plnSamp.multScen.relRangeShift;
        caSampRes(i).absRangeShift = plnSamp.multScen.absRangeShift;
        caSampRes(i).isoShift      = plnSamp.multScen.isoShift;

        if isfield(resultSamp,[plnSamp.bioParam.quantityVis '_deformed'])
            caSampRes(i).dvh = matRad_calcDVH(cstSamp,resultSamp.(plnSamp.bioParam.quantityVis),'cum',dvhPoints);
            caSampRes(i).qi  = matRad_calcQualityIndicators(cstSamp,plnSamp,resultSamp.(plnSamp.bioParam.quantityVis)*plnSamp.numOfFractions,refGy,refVol);
        else
            caSampRes(i).dvh = matRad_calcDVH(cstSamp,resultSamp.([plnSamp.bioParam.quantityVis '_deformed']),'cum',dvhPoints);
            caSampRes(i).qi  = matRad_calcQualityIndicators(cstSamp,plnSamp,resultSamp.([plnSamp.bioParam.quantityVis '_deformed'])*plnSamp.numOfFractions,refGy,refVol);
        end
        
    end
    
end

%% add subindices
pln.subIx        = subIx;

end
