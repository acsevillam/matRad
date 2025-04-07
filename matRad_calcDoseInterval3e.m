function [dij_dummy, pln_dummy, dij, pln, dij_interval] = matRad_calcDoseInterval3e(ct, cst, stf, pln, dij, targetStructSel, OARStructSel, kdin, kmax, retentionThreshold, file_pattern)

matRad_cfg = MatRad_Config.instance();
[env, envver] = matRad_getEnvironment();
pln_dummy = pln;

if ~exist('kdin','var') || isempty(kdin)
    kdin='dinamic';
end

if ~exist('kmax','var') || isempty(kmax)
    kmax=10;
end

if ~exist('retentionThreshold','var') || isempty(retentionThreshold)
    retentionThreshold=1.0;
end

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

% Attempt to load precomputed dij_interval.center and dij_interval.radius
[center_loaded, radius_loaded, file_loaded] = matRad_loadTargetDij(file_pattern, targetSubIx);

if ~isempty(center_loaded) && ~isempty(radius_loaded)
    fprintf('[matRad_calcDoseInterval3e] Loaded precomputed center and radius from: %s\n', file_loaded);
    dij_interval.center(targetSubIx, :) = center_loaded;
    dij_interval.radius = radius_loaded;
else
    % Estimate the memory usage per voxel in MB (empirically estimated)
    memory_factor = 0.05;
    estimatedMemoryPerVoxelMB = (dij.totalNumOfBixels+dij.totalNumOfBixels*dij.totalNumOfBixels)*8/1E6 * memory_factor;
    
    % Get available RAM depending on OS
    availableMB = matRad_getAvailableMemoryMB();
    
    % Compute the memory budget for all batches (e.g., 80% of available memory)
    totalBatchMemoryBudgetMB = 0.80 * availableMB;
    
    % Adjust the maximum batch size depending on memory per voxel
    maxBatchSize = floor(totalBatchMemoryBudgetMB / estimatedMemoryPerVoxelMB / nWorkers);
    
    % Print estimated configuration
    fprintf('Available RAM: %.1f MB | Estimated per-voxel: %.2f MB | nWorkers: %d | maxBatchSize: %d voxels\n', ...
        availableMB, estimatedMemoryPerVoxelMB, nWorkers, maxBatchSize);
    
    nBatches=1;
    limitBatches = 1000;
    
    % Dynamically increase nBatches if the batch size exceeds the memory-safe threshold
    while (ceil(numel(targetSubIx) / nBatches) > maxBatchSize) || nBatches == limitBatches
        nBatches = nBatches + 1;
    end
    
    % Recalculate target batch size after adjustment
    targetBatchSize = ceil(numel(targetSubIx) / nBatches);
    
    % Print and log results
    logMsg = sprintf(['[matRad batching log]\nTimestamp: %s\n', 'Available RAM: %.1f MB\nEstimated per-voxel: %.2f MB\n',...
        'nWorkers: %d\nEstimated maxBatchSize: %d voxels\n', 'Final nBatches: %d | targetBatchSize: %d voxels\n\n'], ...
        datestr(now), availableMB, estimatedMemoryPerVoxelMB, nWorkers, maxBatchSize, nBatches, targetBatchSize);
    
    fprintf(logMsg);
    
    if exist('parfor_progress.txt', 'file') ~= 2
        fclose(fopen('parfor_progress.txt', 'w'));
    end
    
    for b = 1:nBatches
        idx_start = (b-1)*targetBatchSize + 1;
        idx_end = min(b*targetBatchSize, numel(targetSubIx));
        currentBatch = targetSubIx(idx_start:idx_end);
    
        fprintf('Processing batch %d of %d (%d voxels)\n', b, nBatches, numel(currentBatch));
        if exist('parfor_progress', 'file') == 2
            FlagParforProgressDisp = true;
            parfor_progress(round(numel(currentBatch)/10));
        else
            FlagParforProgressDisp = false;
        end
    
        numBixels = size(dij_list{1}, 2);
        numVoxelsInBatch = numel(currentBatch);
    
        % Reduce dose list for the current batch
        dij_list_reduced = cellfun(@(data) data(currentBatch,:), dij_list, 'UniformOutput', false);
        numScenarios = numel(dij_list_reduced);
    
        % Preallocate result matrices
        centers_block = zeros(numVoxelsInBatch, numBixels);
        radius_block_local = zeros(numBixels, numBixels);  % Accumulated radius
    
        % Vectorized parfor
        parfor it = 1:numVoxelsInBatch
            dij_tmp = zeros(numScenarios, numBixels);
    
            for s = 1:numScenarios
                dij_tmp(s, :) = dij_list_reduced{s}(it, :);
            end
    
            % Weighted contributions
            dij_tmp_weighted = dij_tmp .* scenProb;
    
            % Center (E[d])
            center_voxel = sum(dij_tmp_weighted, 1);
    
            % Radius (E[d^2])
            radius_voxel = dij_tmp' * dij_tmp_weighted;
    
            % Store result
            centers_block(it, :) = center_voxel;
            
            % Accumulate radius using reduction strategy
            radius_block_local = radius_block_local + radius_voxel;
    
            % Optional progress
            if FlagParforProgressDisp && mod(it,10)==0
                parfor_progress;
            end
        end
    
        if FlagParforProgressDisp
            parfor_progress(0);
        end

        fprintf('Finishing batch calculation ... \n');
        toc
    
        % Assign result
        dij_interval.center(currentBatch, :) = centers_block;
        dij_interval.radius = dij_interval.radius + radius_block_local;
    
        clear dij_list_reduced;
        fprintf('Finishing batch %d storage ...\n', b);
        toc
    end

end

whos dij_interval;
toc

% Inicia el parpool solo si no está abierto
if isempty(gcp('nocreate'))
    parpool('local', nWorkers);
end

dij_interval.OARSubIx = OARSubIx;

% OAR voxel batching
tic

% Estimate the memory usage per OAR voxel (can be slightly higher due to SVD storage)
memory_factor = 1;
estimatedMemoryPerOARVoxelMB = (kmax*kmax + 2*kmax*dij.totalNumOfBixels)*8/1E6 * memory_factor;

% Get available RAM depending on OS
availableMB = matRad_getAvailableMemoryMB();

% Compute memory budget per worker (e.g., 80% of available memory)
totalOARBatchMemoryBudgetMB = 0.80 * availableMB;
maxOARBatchSize = floor(totalOARBatchMemoryBudgetMB / estimatedMemoryPerOARVoxelMB / nWorkers);

% Initial estimate based on available workers
nOARBatches = 1;
limitOARBatches = 100;

% Dynamically increase nOARBatches if the batch size exceeds the memory-safe threshold
while (ceil(numel(OARSubIx) / nOARBatches) > maxOARBatchSize) || nOARBatches == limitOARBatches
    nOARBatches = nOARBatches + 1;
end

% Recalculate batch size after adjustment
OARBatchSize = ceil(numel(OARSubIx) / nOARBatches);

% Print and log results
logMsg = sprintf(['[matRad OAR batching log]\nTimestamp: %s\n', 'Available RAM: %.1f MB\nEstimated per-OAR-voxel: %.2f MB\n',...
    'nWorkers: %d\nEstimated maxOARBatchSize: %d voxels\n','Final nOARBatches: %d | OARBatchSize: %d voxels\n\n'], ...
    datestr(now), availableMB, estimatedMemoryPerOARVoxelMB, nWorkers, maxOARBatchSize, nOARBatches, OARBatchSize);

fprintf(logMsg);

if exist('parfor_progress.txt', 'file') ~= 2
    fclose(fopen('parfor_progress.txt', 'w'));
end

for b = 1:nOARBatches
    idx_start = (b-1)*OARBatchSize + 1;
    idx_end = min(b*OARBatchSize, numel(OARSubIx));
    currentBatch = OARSubIx(idx_start:idx_end);
    
    fprintf('Processing batch %d of %d (%d voxeles) \n', b, nOARBatches, numel(currentBatch));
    if exist('parfor_progress', 'file') == 2
        FlagParforProgressDisp = true;
        parfor_progress(round(numel(currentBatch)/10));  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
    else
        matRad_cfg.dispInfo('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n');
        FlagParforProgressDisp = false;
    end

    % Outside parfor
    numBixels = size(dij_list{1}, 2);
    numOARVoxelsInBatch = numel(currentBatch);

    dij_list_reduced = cellfun(@(data) data(currentBatch,:), dij_list, 'UniformOutput', false);
    dij_batch_OAR = repmat(struct('Ix', [], 'center', [], 'k', [], 'U', [], 'S', [], 'V', []), numel(currentBatch), 1);

    parfor it = 1:numel(currentBatch)
    
        % Initialize temporary matrix to store dose values across scenarios for this voxel
        dij_tmp = zeros(numel(dij_list_reduced), size(dij_list_reduced{1}, 2));  % (numScenarios x numBixels)
    
        % Fill dij_tmp with the dose contributions for voxel 'it' across all scenarios
        for s = 1:numel(dij_list_reduced)
            dij_tmp(s, :) = dij_list_reduced{s}(it, :);  % Row 'it' from scenario 's'
        end
    
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
           
        % Select k according to kdin
        if isequal(kdin,'dinamic')
            % Perform truncated Singular Value Decomposition (SVD)
            [U, S, V] = svds(covMatrix, kmax, 'largest');  % Keep top 10 singular values/vectors
        
            % Extract singular values and compute total and cumulative energy
            singularValues = diag(S);
            totalEnergy = sum(singularValues.^2);
            cumulativeEnergy = cumsum(singularValues.^2);
        
            % Find number of components 'k' to retain enough energy (based on threshold)
            k = find(cumulativeEnergy / totalEnergy >= retentionThreshold, 1, 'first');

            % Store reduced-rank SVD components in sparse format
            dij_batch_OAR(it).k = k;
            dij_batch_OAR(it).U = sparse(U(:, 1:k));
            dij_batch_OAR(it).S = sparse(S(1:k, 1:k));
            dij_batch_OAR(it).V = sparse(V(:, 1:k));

        elseif isequal(kdin,'static')
            k=kmax;
            % Perform truncated Singular Value Decomposition (SVD)
            [U, S, V] = svds(covMatrix, k, 'largest');
            % Store reduced-rank SVD components in sparse format
            dij_batch_OAR(it).k = k;
            dij_batch_OAR(it).U = sparse(U);
            dij_batch_OAR(it).S = sparse(S);
            dij_batch_OAR(it).V = sparse(V);            
        else
            disp('k dinamics did not find!');
        end
    
        % Optional progress display inside the parfor loop
        if FlagParforProgressDisp && mod(it,10)==0
            %fprintf('Processing batch %d of %d \n', b, nOARBatches);
            parfor_progress;
        end
    end

    if FlagParforProgressDisp
        parfor_progress(0);
    end

    % Free memory (optional but useful in large cases)
    clear dij_list_reduced;

    fprintf('Finishing batch calculation ... \n');
    toc
    % Vectorized accumulation at the end of each batch

    % Preallocate center and radius blocks for batch
    centers_block = zeros(numOARVoxelsInBatch, numBixels);      % (voxels x bixels)
    k_block = zeros(numOARVoxelsInBatch, 1);      % (voxels x bixels)
    U_block = cell(numOARVoxelsInBatch, 1);
    S_block = cell(numOARVoxelsInBatch, 1);
    V_block = cell(numOARVoxelsInBatch, 1);

    for it = 1:numOARVoxelsInBatch
        centers_block(it, :) = dij_batch_OAR(it).center;         % Store each voxel center
        k_block(it, :) = dij_batch_OAR(it).k;
        U_block{it} = dij_batch_OAR(it).U;
        S_block{it} = dij_batch_OAR(it).S;
        V_block{it} = dij_batch_OAR(it).V;
    end

    % Assign in block
    dij_interval.center(currentBatch, :) = centers_block;
    dij_interval.k(currentBatch, :) = k_block;
    dij_interval.U(idx_start:idx_end) = U_block;
    dij_interval.S(idx_start:idx_end) = S_block;
    dij_interval.V(idx_start:idx_end) = V_block;
    
    clear dij_batch_OAR;
    fprintf('Finishing batch data storage ... \n');
    toc

end

whos dij_interval;
toc

end

function [center, radius, fileName] = matRad_loadTargetDij(file_pattern, Ix)
% matRad_loadTargetDij - Loads dij_interval.center and radius from the first matching file.
%
% Inputs:
%   file_pattern : wildcard string pattern to match (e.g., 'folder/*.mat')
%   Ix           : voxel indices to extract from dij_interval.center
%
% Outputs:
%   center       : rows of dij_interval.center corresponding to Ix
%   radius       : dij_interval.radius (full matrix)
%   fileName     : full path of the loaded file, empty if none was loaded

    center = [];
    radius = [];
    fileName = '';
    
    files = dir(file_pattern);
    
    if isempty(files)
        fprintf('[matRad_loadTargetDij] No file found matching pattern: %s\n', file_pattern);
        return;
    end

    % Take the first match
    fileName = fullfile(files(1).folder, files(1).name);
    fprintf('[matRad_loadTargetDij] Loading file: %s\n', fileName);
    
    try
        f = load(fileName, 'dij_interval');
        
        if isfield(f, 'dij_interval')
            dij = f.dij_interval;

            % Extract center if available and Ix is valid
            if isfield(dij, 'center') && ~isempty(Ix)
                center = dij.center(Ix, :);
            else
                warning('[matRad_loadTargetDij] Field ''dij_interval.center'' missing or Ix invalid in %s.', fileName);
                center = [];
                fileName = '';
            end

            % Extract radius if available
            if isfield(dij, 'radius')
                radius = dij.radius;
            else
                warning('[matRad_loadTargetDij] Field ''dij_interval.radius'' not found in %s.', fileName);
                radius = [];
                fileName = '';
            end
        else
            warning('[matRad_loadTargetDij] ''dij_interval'' structure not found in %s.', fileName);
            center = [];
            radius = [];
            fileName = '';
        end

    catch ME
        warning('[matRad_loadTargetDij] Failed to load file: %s\nError: %s', fileName, ME.message);
        center = [];
        radius = [];
        fileName = '';
    end
end

function availableMB = matRad_getAvailableMemoryMB()
% matRad_getAvailableMemoryMB - Returns available system memory in MB
% Works on Windows, macOS, and Linux (including Rocky).
%
% Output:
%   availableMB : Approximate available RAM in MB.

    if ispc
        user = memory;
        availableMB = user.MemAvailableAllArrays / 1e6;

    elseif ismac
        [~, memInfo] = system('vm_stat');
        pageSizeLine = regexp(memInfo, 'page size of (\d+) bytes', 'tokens', 'once');
        if isempty(pageSizeLine)
            warning('Could not parse vm_stat output.');
            availableMB = 0;
            return;
        end
        pageSize = str2double(pageSizeLine{1});
        freePages = sum(cellfun(@(x) sscanf(x{2}, '%d'), ...
            regexp(memInfo, 'Pages (free|inactive|speculative):\s+(\d+)', 'tokens')));
        availableMB = (freePages * pageSize) / 1e6;

    elseif isunix
        [~, memInfo] = system('grep MemAvailable /proc/meminfo');
        tokens = regexp(memInfo, 'MemAvailable:\s+(\d+)\s+kB', 'tokens');
        if ~isempty(tokens)
            availableMB = str2double(tokens{1}) / 1024;
        else
            warning('Could not read /proc/meminfo.');
            availableMB = 0;
        end
    else
        error('Unsupported operating system.');
    end
end
