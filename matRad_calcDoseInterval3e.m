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

nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));

% Fallback para pruebas locales
if isnan(nWorkers) || nWorkers == 0
    nCores = feature('numcores');
    nWorkers = max(1, nCores - 2);
end

% Inicia el parpool solo si no está abierto
if isempty(gcp('nocreate'))
    parpool('local', nWorkers);
end

dij_interval.center = sparse(dij.doseGrid.numOfVoxels, dij.totalNumOfBixels);
dij_interval.radius = sparse(dij.totalNumOfBixels, dij.totalNumOfBixels);
dij_interval.targetSubIx = targetSubIx;

% Target voxel batching
tic
nBatches = min([ceil(numel(targetSubIx)/nWorkers/2) 100]);
targetBatchSize = ceil(numel(targetSubIx) / nBatches);

if exist('parfor_progress.txt', 'file') ~= 2
    fclose(fopen('parfor_progress.txt', 'w'));
end

for b = 1:nBatches
    
    idx_start = (b-1)*targetBatchSize + 1;
    idx_end = min(b*targetBatchSize, numel(targetSubIx));
    currentBatch = targetSubIx(idx_start:idx_end);

    dij_batch = repmat(struct('Ix', [], 'center', [], 'radius', []), numel(currentBatch), 1);

    fprintf('Processing batch %d of %d (%d voxeles)', b, nBatches, numel(currentBatch));
    if exist('parfor_progress', 'file') == 2
        FlagParforProgressDisp = true;
        parfor_progress(round(numel(currentBatch)/10));  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
    else
        matRad_cfg.dispInfo('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n');
        FlagParforProgressDisp = false;
    end

    dij_list_reduced = cellfun(@(data) data(currentBatch,:), dij_list, 'UniformOutput', false);

    parfor it = 1:numel(currentBatch)
        Ix = currentBatch(it); % Global voxel index being processed
    
        % Extract the dose values across all scenarios for this voxel
        dij_tmp = cell2mat(cellfun(@(data) data(it,:), dij_list_reduced, 'UniformOutput', false));  
        % dij_tmp: (numScenarios x numBixels)
    
        % Apply scenario probabilities to the dose values (weighted)
        dij_tmp_weighted = dij_tmp .* scenProb;  % (numScenarios x numBixels)
    
        % Store voxel index
        dij_batch(it).Ix = Ix;
    
        % Compute the expected dose per bixel for this voxel (center of interval)
        dij_batch(it).center = sum(dij_tmp_weighted, 1)';  % (numBixels x 1)
    
        % Compute E[d^2] per bixel for this voxel
        dij_batch(it).radius = (dij_tmp' * dij_tmp_weighted)';  % (numBixels x 1)
        % Note: dij_tmp' * dij_tmp_weighted = sum_k w_k * d_k^2 (per bixel)
    
        % Optional progress display
        if FlagParforProgressDisp && mod(it,10)==0
            fprintf('Processing batch %d of %d', b, nBatches);
            parfor_progress;
        end
    end
    
    % Free memory (optional but useful in large cases)
    clear dij_list_reduced dij_tmp dij_tmp_weighted;
    
    % Accumulate results from each voxel
    for it = 1:numel(currentBatch)
        % Store expected dose (center) in output matrix
        dij_interval.center(currentBatch(it), :) = dij_batch(it).center;
    
        % Accumulate the E[d^2] terms across voxels (to later compute total variance)
        dij_interval.radius = dij_interval.radius + dij_batch(it).radius;
    end

    if FlagParforProgressDisp
        parfor_progress(0);
    end

    clear dij_batch;
    toc
end

whos dij_interval;
toc

% OAR voxel batching
tic
dij_interval.OARSubIx = OARSubIx;
nOARBatches = min([ceil(numel(OARSubIx)/nWorkers/2) 100]);
OARBatchSize = ceil(numel(OARSubIx) / nOARBatches);

if exist('parfor_progress.txt', 'file') ~= 2
    fclose(fopen('parfor_progress.txt', 'w'));
end

for b = 1:nOARBatches
    idx_start = (b-1)*OARBatchSize + 1;
    idx_end = min(b*OARBatchSize, numel(OARSubIx));
    currentBatch = OARSubIx(idx_start:idx_end);
    
    fprintf('Processing batch %d of %d (%d voxeles)', b, nOARBatches, numel(currentBatch));
    if exist('parfor_progress', 'file') == 2
        FlagParforProgressDisp = true;
        parfor_progress(round(numel(currentBatch)/10));  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
    else
        matRad_cfg.dispInfo('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n');
        FlagParforProgressDisp = false;
    end

    dij_list_reduced = cellfun(@(data) data(currentBatch,:), dij_list, 'UniformOutput', false);
    dij_batch_OAR = repmat(struct('Ix', [], 'center', [], 'U', [], 'S', [], 'V', []), numel(currentBatch), 1);

    parfor it = 1:numel(currentBatch)
        Ix = currentBatch(it);         % Global index of the voxel being processed
    
        % Initialize temporary matrix to store dose values across scenarios for this voxel
        dij_tmp = zeros(numel(dij_list_reduced), size(dij_list_reduced{1}, 2));  % (numScenarios x numBixels)
    
        % Fill dij_tmp with the dose contributions for voxel 'it' across all scenarios
        for s = 1:numel(dij_list_reduced)
            dij_tmp(s, :) = dij_list_reduced{s}(it, :);  % Row 'it' from scenario 's'
        end
    
        % Store voxel index in output structure
        dij_batch_OAR(it).Ix = Ix;
    
        % Compute the expected dose (center of the dose interval) across scenarios
        dij_batch_OAR(it).center = sum(dij_tmp .* scenProb, 1);  % (1 x numBixels)
    
        % Center the dose matrix by subtracting the expected dose
        dij_centered = dij_tmp - dij_batch_OAR(it).center;  % (numScenarios x numBixels)
    
        % Apply square root of scenario probabilities (for weighted covariance)
        dij_centered_weighted = dij_centered .* sqrt(scenProb);  % (numScenarios x numBixels)
    
        % Compute the weighted covariance matrix for the bixels
        covMatrix = dij_centered_weighted' * dij_centered_weighted;  % (numBixels x numBixels)
        %covMatrix=(dij_tmp'*diag(scenProb)*dij_tmp-dij_batch_OAR(it).center'*dij_batch_OAR(it).center); % (numBixels x numBixels)
        % Note: covMatrix = E[d^2]-E[d]^2
    
        % Perform truncated Singular Value Decomposition (SVD)
        [U, S, V] = svds(covMatrix, 10, 'largest');  % Keep top 10 singular values/vectors
    
        % Extract singular values and compute total and cumulative energy
        singularValues = diag(S);
        totalEnergy = sum(singularValues.^2);
        cumulativeEnergy = cumsum(singularValues.^2);
    
        % Find number of components 'k' to retain enough energy (based on threshold)
        k = find(cumulativeEnergy / totalEnergy >= retentionThreshold, 1, 'first');
    
        % Store reduced-rank SVD components in sparse format
        dij_batch_OAR(it).U = sparse(U(:, 1:k));
        dij_batch_OAR(it).S = sparse(S(1:k, 1:k));
        dij_batch_OAR(it).V = sparse(V(:, 1:k));
    
        % Optional progress display inside the parfor loop
        if FlagParforProgressDisp && mod(it,10)==0
            fprintf('Processing batch %d of %d', b, nOARBatches);
            parfor_progress;
        end
    end

    if FlagParforProgressDisp
        parfor_progress(0);
    end

    % Accumulate results from each voxel
    for it=1:numel(currentBatch)
        % Store expected dose (center) in output matrix
        dij_interval.center(currentBatch(it),:)=dij_batch_OAR(it).center;

        % Store [U,S,V] for covMatrix = E[d^2]-E[d]^2 in output matrix
        dij_interval.U{idx_start+it-1}=dij_batch_OAR(it).U;
        dij_interval.S{idx_start+it-1}=dij_batch_OAR(it).S;
        dij_interval.V{idx_start+it-1}=dij_batch_OAR(it).V;
    end

    clear dij_batch_OAR;
    toc
end

whos dij_interval;
toc

end
