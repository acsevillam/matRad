function [dij_dummy, pln_dummy, dij, pln, dij_interval] = matRad_calcDoseInterval3d(ct, cst, stf, pln, dij, targetStructSel, OARStructSel, retentionThreshold)

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
clear cst;

center_local = cell(1, numel(targetSubIx));
radius_local = cell(1, numel(targetSubIx));

parfor it = 1:numel(targetSubIx)
    Ix = targetSubIx(it);
    dij_tmp = cell2mat(cellfun(@(data) data(Ix,:), dij_list, 'UniformOutput', false));
    dij_tmp_weighted = dij_tmp .* scenProb';
    center_local{it} = sum(dij_tmp_weighted, 1)';
    radius_local{it} = dij_tmp' * dij_tmp_weighted;
end
clear scenIx;

% Assemble dij_interval
nBixels = dij.totalNumOfBixels;
dij_interval.center = sparse(dij.doseGrid.numOfVoxels, nBixels);
dij_interval.radius = sparse(nBixels, nBixels);
dij_interval.targetSubIx = targetSubIx;

for it = 1:numel(targetSubIx)
    dij_interval.center(targetSubIx(it), :) = center_local{it};
    dij_interval.radius = dij_interval.radius + radius_local{it};
end
clear center_local;
clear radius_local;

batch_size = 100;

% Define batch indices
nOAR = numel(OARSubIx);
batch_indices = cell(ceil(nOAR / batch_size), 1);
OAR_batches = cell(ceil(nOAR / batch_size), 1);
for b = 1:numel(batch_indices)
    batch_indices{b} = ((b-1)*batch_size + 1):min(b*batch_size, nOAR);
    OAR_batches{b} = OARSubIx(batch_indices{b});
end

% Initialize containers
dij_interval.U = cell(size(OARSubIx));
dij_interval.S = cell(size(OARSubIx));
dij_interval.V = cell(size(OARSubIx));
dij_interval.OARSubIx = OARSubIx;

batch_results = cell(numel(batch_indices), 1);

parfor b = 1:numel(batch_indices)
    OAR_batch = OAR_batches{b};

    center_oar_local = cell(1, numel(OAR_batch));
    U_local = cell(1, numel(OAR_batch));
    S_local = cell(1, numel(OAR_batch));
    V_local = cell(1, numel(OAR_batch));

    for idx = 1:numel(OAR_batch)
        Ix = OAR_batch(idx);
        dij_tmp = cell2mat(cellfun(@(data) data(Ix,:), dij_list, 'UniformOutput', false));
        dij_tmp_weighted = dij_tmp .* scenProb';
        center_tmp = sum(dij_tmp_weighted, 1)';
        radius_tmp = sparse(dij_tmp' * dij_tmp_weighted - center_tmp * center_tmp');

        [U, S, V] = svds(radius_tmp, 10, 'largest');
        singularValues = diag(S);
        totalEnergy = sum(singularValues.^2);
        cumulativeEnergy = cumsum(singularValues.^2);
        k = find(cumulativeEnergy / totalEnergy >= retentionThreshold, 1, 'first');

        center_oar_local{idx} = center_tmp;
        U_local{idx} = sparse(U(:, 1:k));
        S_local{idx} = sparse(S(1:k, 1:k));
        V_local{idx} = sparse(V(:, 1:k));
    end

    batch_results{b} = struct('it', batch_indices{b}, 'center', {center_oar_local}, ...
                              'U', {U_local}, 'S', {S_local}, 'V', {V_local});
end

for b = 1:numel(batch_results)
    it_batch = batch_results{b}.it;
    for idx = 1:numel(it_batch)
        it = it_batch(idx);
        dij_interval.center(OARSubIx(it), :) = batch_results{b}.center{idx};
        dij_interval.U{it} = batch_results{b}.U{idx};
        dij_interval.S{it} = batch_results{b}.S{idx};
        dij_interval.V{it} = batch_results{b}.V{idx};
    end
end
clear dij_list;

end
