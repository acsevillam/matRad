function [dij_dummy, pln_dummy, dij, pln, dij_interval] = matRad_calcDoseInterval3c(ct, cst, stf, pln, dij, targetStructSel, OARStructSel, retentionThreshold)

matRad_cfg = MatRad_Config.instance();
[env, envver] = matRad_getEnvironment();
pln_dummy = pln;

% Retrieve 2 dummy case scenarios for dose calculation and optimization
pln_dummy.multScen = matRad_multScen(ct, 'nomScen');

% Apply deformation vector field if available
if isfield(ct, 'dvf')
    if ~isequal(dij.doseGrid.dimensions, ct.dvfDim) || ~isequal(ct.dvfType, 'push')
        matRad_cfg.dispWarning('Dose cube and deformation vector field dimensions are not equal. \n');
        
        % Instantiate elastic registration
        metadata.nItera = 100;
        metadata.dvfType = 'push';
        register = matRad_ElasticImageRegistration(ct, cst, 1, metadata);
        clear 'metadata';
        
        % Resize dvf according to dose grid
        [ct] = register.calcDVF_resized(dij.doseGrid.resolution);
    end
end

% Calculate dummy case dij to save interval
dij_dummy = matRad_calcPhotonDose(ct, stf, pln_dummy, cst);

scenIx = find(pln.multScen.scenMask);
dij_list = dij.physicalDose(scenIx);
scenProb = pln.multScen.scenProb;

% Resizing cst to dose cube resolution
cst = matRad_resizeCstToGrid(cst, dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);

targetV = [];
% define voxels for sampling
if ~exist('targetStructSel', 'var') || sum(size(targetStructSel)) == 0
    targetV = [cst{:,4}];
    % final voxel subset for sampling
    targetSubIx = unique(vertcat(targetV{:}));
else
    for i=1:size(cst,1)
        for j = 1:numel(targetStructSel)
            if strcmp(targetStructSel{j}, cst{i,2})
                targetV = [targetV; cst{i,4}{1}];
            end
        end
    end
    % final voxel subset for sampling
    targetSubIx = targetV;
end

OARV = [];
% define voxels for sampling
if ~exist('OARStructSel', 'var') || sum(size(OARStructSel)) == 0
    OARV = [cst{:,4}];
    % final voxel subset for sampling
    OARSubIx = unique(vertcat(OARV{:}));
else
    for i=1:size(cst,1)
        for j = 1:numel(OARStructSel)
            if strcmp(OARStructSel{j}, cst{i,2})
                OARV = [OARV; cst{i,4}{1}];
            end
        end
    end
    % final voxel subset for sampling
    OARSubIx = OARV;
end

% Parallel pool check
p = gcp(); % If no pool, create new one.
    
if isempty(p)
    poolSize = 1;
else
    poolSize = p.NumWorkers;
end

dij_interval.center = sparse(dij.doseGrid.numOfVoxels, dij.totalNumOfBixels);
dij_interval.radius = sparse(dij.totalNumOfBixels, dij.totalNumOfBixels);
dij_interval.targetSubIx = targetSubIx;

for it = 1:numel(targetSubIx)
    Ix = targetSubIx(it);
    dij_tmp = cell2mat(cellfun(@(data) data(Ix,:), dij_list, 'UniformOutput', false));
    dij_interval.center(targetSubIx(it), :) = sum(dij_tmp' * diag(scenProb), 2);
    dij_interval.radius = dij_interval.radius + dij_tmp' * diag(scenProb) * dij_tmp;
    if mod(it, 10) == 0 || it == numel(targetSubIx)
        fprintf('Progress: %.2f%%\n', (it / numel(targetSubIx)) * 100);
    end
end


% SVD Decomposition with dynamic k selection

dij_interval.U = cell(size(OARSubIx));
dij_interval.S = cell(size(OARSubIx));
dij_interval.V = cell(size(OARSubIx));
dij_interval.OARSubIx = OARSubIx;

for it = 1:numel(OARSubIx)
    Ix = OARSubIx(it);
    dij_tmp = cell2mat(cellfun(@(data) data(Ix,:), dij_list, 'UniformOutput', false));
    dij_interval.center(OARSubIx(it), :) = sum(dij_tmp' * diag(scenProb), 2);

    radius_tmp = sparse(dij_tmp' * diag(scenProb) * dij_tmp - dij_interval.center(OARSubIx(it), :)' * dij_interval.center(OARSubIx(it), :));
    
    % Compute full SVD
    %[U, S, V] = svd(radius_tmp, 'econ');
    % Compute 5 largest singular elements
    [U, S, V] = svds(radius_tmp,20,'largest'); 
    
    % Determine dynamic k based on retention threshold
    singularValues = diag(S);
    totalEnergy = sum(singularValues.^2);
    cumulativeEnergy = cumsum(singularValues.^2);
    k = find(cumulativeEnergy / totalEnergy >= retentionThreshold, 1, 'first');
    
    % Store only the necessary components
    dij_interval.U{it} = sparse(U(:, 1:k));
    dij_interval.S{it} = sparse(S(1:k, 1:k));
    dij_interval.V{it} = sparse(V(:, 1:k));

    if mod(it, 10) == 0 || it == numel(OARSubIx)
        fprintf('Progress: %.2f%%\n', (it / numel(OARSubIx)) * 100);
    end
end

end
