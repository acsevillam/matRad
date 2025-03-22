function [dij_dummy, pln_dummy, dij, pln, dij_interval] = matRad_calcDoseInterval3e(ct, cst, stf, pln, dij, targetStructSel, OARStructSel, retentionThreshold)

matRad_cfg = MatRad_Config.instance();
[env, envver] = matRad_getEnvironment();
pln_dummy = pln;

% Retrieve dummy case scenario
dummyScen = matRad_multScen(ct, 'nomScen');
pln_dummy.multScen = dummyScen;

% Apply deformation vector field if available
if isfield(ct, 'dvf')
    if ~isequal(dij.doseGrid.dimensions, ct.dvfDim) || ~isequal(ct.dvfType, 'push')
        matRad_cfg.dispWarning('Dose cube and deformation vector field dimensions are not equal. \n');
        metadata.nItera = 100;
        metadata.dvfType = 'push';
        register = matRad_ElasticImageRegistration(ct, cst, 1, metadata);
        clear metadata;
        [ct] = register.calcDVF_resized(dij.doseGrid.resolution);
        clear register;
    end
end

% Calculate dummy dose
dij_dummy = matRad_calcPhotonDose(ct, stf, pln_dummy, cst);

% Scenario probabilities and dose list
scenIx = find(pln.multScen.scenMask);
dij_list = dij.physicalDose(scenIx);
scenProb = pln.multScen.scenProb;

% Resize cst to dose grid
cst = matRad_resizeCstToGrid(cst, dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);

% Target voxel selection
if ~exist('targetStructSel', 'var') || sum(size(targetStructSel)) == 0
    targetV = [cst{:,4}];
    targetSubIx = unique(vertcat(targetV{:}));
else
    targetV = cell(size(cst,1) * numel(targetStructSel), 1);
    counter = 0;
    for i = 1:size(cst,1)
        for j = 1:numel(targetStructSel)
            if strcmp(targetStructSel{j}, cst{i,2})
                counter = counter + 1;
                targetV{counter} = cst{i,4}{1};
            end
        end
    end
    targetV = vertcat(targetV{1:counter});
    targetSubIx = targetV;
end
clear targetV;
clear counter;

% OAR voxel selection
if ~exist('OARStructSel', 'var') || sum(size(OARStructSel)) == 0
    OARV = [cst{:,4}];
    OARSubIx = unique(vertcat(OARV{:}));
else
    OARV = cell(size(cst,1) * numel(OARStructSel), 1);
    counter = 0;
    for i = 1:size(cst,1)
        for j = 1:numel(OARStructSel)
            if strcmp(OARStructSel{j}, cst{i,2})
                counter = counter + 1;
                OARV{counter} = cst{i,4}{1};
            end
        end
    end
    OARV = vertcat(OARV{1:counter});
    OARSubIx = OARV;
end
clear OARV;
clear counter;
clear cst;

nCores = feature('numcores');
nWorkers = max(1, nCores - 2); % deja libres 2 núcleos

if isempty(gcp('nocreate'))
    parpool('local', nWorkers);
end

dij_interval.center = sparse(dij.doseGrid.numOfVoxels, dij.totalNumOfBixels);
dij_interval.radius = sparse(dij.totalNumOfBixels, dij.totalNumOfBixels);
dij_interval.targetSubIx = targetSubIx;

% Target voxel batching
nBatches = 100;
batch_size = ceil(numel(targetSubIx) / nBatches);
for b = 1:nBatches
    idx_start = (b-1)*batch_size + 1;
    idx_end = min(b*batch_size, numel(targetSubIx));
    currentBatch = targetSubIx(idx_start:idx_end);

    dij_batch = repmat(struct('Ix', [], 'center', [], 'radius', []), numel(currentBatch), 1);

    if exist('parfor_progress', 'file') == 2
        FlagParforProgressDisp = true;
        parfor_progress(round(numel(targetSubIx(idx_start:idx_end))/10000));  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
    else
        matRad_cfg.dispInfo('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n');
        FlagParforProgressDisp = false;
    end

    parfor it = 1:numel(currentBatch)
        Ix = currentBatch(it);
        dij_tmp = cell2mat(cellfun(@(data) data(Ix,:), dij_list, 'UniformOutput', false));
        dij_tmp_weighted = dij_tmp .* scenProb;
        dij_batch(it).Ix = Ix;
        dij_batch(it).center = sum(dij_tmp_weighted, 1);
        dij_batch(it).radius = dij_tmp' * dij_tmp_weighted;

        if FlagParforProgressDisp && mod(it,10000)==0
            parfor_progress;
        end
    end

    for it = 1:numel(currentBatch)
        dij_interval.center(currentBatch(it), :) = dij_batch(it).center;
        dij_interval.radius = dij_interval.radius + dij_batch(it).radius;
    end

    if FlagParforProgressDisp
        parfor_progress(0);
    end

    clear dij_batch;
end

% Preallocate structure
dij_interval_OAR = repmat(struct('Ix', [], 'center', [], 'U', [], 'S', [], 'V', []), numel(OARSubIx), 1);

if exist('parfor_progress.txt', 'file') ~= 2
    fclose(fopen('parfor_progress.txt', 'w'));
end

if exist('parfor_progress', 'file') == 2
    FlagParforProgressDisp = true;
    parfor_progress(round(numel(OARSubIx)/10000));  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
else
    matRad_cfg.dispInfo('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n');
    FlagParforProgressDisp = false;
end

parfor it=1:numel(OARSubIx)

    Ix=OARSubIx(it);
    dij_tmp=cell2mat(cellfun(@(data) data(Ix,:),dij_list,'UniformOutput',false));
    dij_tmp_weighted = dij_tmp .* scenProb; % vector-matrix element-wise

    % Voxel index
    dij_interval_OAR(it).Ix=Ix;

    % Interval center dose influence matrix
    dij_interval_OAR(it).center=sum(dij_tmp_weighted, 1);
    
    % Interval radius dose influence matrix
    radius_tmp=(dij_tmp' * dij_tmp_weighted - dij_interval_OAR(it).center * dij_interval_OAR(it).center');


    [U, S, V] = svds(radius_tmp, 10, 'largest');
    singularValues = diag(S);
    totalEnergy = sum(singularValues.^2);
    cumulativeEnergy = cumsum(singularValues.^2);
    k = find(cumulativeEnergy / totalEnergy >= retentionThreshold, 1, 'first');

    dij_interval_OAR(it).U = sparse(U(:, 1:k));
    dij_interval_OAR(it).S = sparse(S(1:k, 1:k));
    dij_interval_OAR(it).V = sparse(V(:, 1:k));

    if FlagParforProgressDisp && mod(it,10000)==0
        parfor_progress;
    end

end

if FlagParforProgressDisp
    parfor_progress(0);
end

dij_interval.U=cell(size(OARSubIx));
dij_interval.S=cell(size(OARSubIx));
dij_interval.V=cell(size(OARSubIx));
dij_interval.OARSubIx=OARSubIx;

for it=1:numel(OARSubIx)
    % Interval center dose influence matrix
    dij_interval.center(OARSubIx(it),:)=dij_interval_OAR(it).center;
    % Interval radius dose influence matrix
    dij_interval.U{it}=dij_interval_OAR(it).U;
    dij_interval.S{it}=dij_interval_OAR(it).S;
    dij_interval.V{it}=dij_interval_OAR(it).V;
end

clear 'dij_interval_OAR';

end
