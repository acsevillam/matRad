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
nWorkers = max(1, nCores - 2);

if isempty(gcp('nocreate'))
    parpool('local', nWorkers);
end

center_local = cell(1, numel(targetSubIx));
radius_local = cell(1, numel(targetSubIx));

parfor it = 1:numel(targetSubIx)
    Ix = targetSubIx(it);
    dij_tmp = cell2mat(cellfun(@(data) data(Ix,:), dij_list, 'UniformOutput', false));
    dij_tmp_weighted = dij_tmp .* scenProb';
    center_tmp = sum(dij_tmp_weighted, 1)';
    radius_tmp = dij_tmp' * dij_tmp_weighted;
    center_local{it} = center_tmp;
    radius_local{it} = radius_tmp;
    
end

for it = 1:numel(targetSubIx)
    dij_interval.center(targetSubIx(it), :) = center_local{it};
    dij_interval.radius = dij_interval.radius + radius_local{it};
end
clear center_local;
clear radius_local;

center_oar = cell(1, numel(OARSubIx));
U_oar = cell(1, numel(OARSubIx));
S_oar = cell(1, numel(OARSubIx));
V_oar = cell(1, numel(OARSubIx));

parfor it = 1:numel(OARSubIx)
    Ix = OARSubIx(it);
    dij_tmp = cell2mat(cellfun(@(data) data(Ix,:), dij_list, 'UniformOutput', false));
    dij_tmp = gpuArray(dij_tmp);
    dij_tmp_weighted = dij_tmp .* gpuArray(scenProb');
    center_tmp = gather(sum(dij_tmp_weighted, 1)');
    radius_tmp = sparse(gather(dij_tmp' * dij_tmp_weighted - center_tmp * center_tmp'));

    [U, S, V] = svds(radius_tmp, 10, 'largest');
    singularValues = diag(S);
    totalEnergy = sum(singularValues.^2);
    cumulativeEnergy = cumsum(singularValues.^2);
    k = find(cumulativeEnergy / totalEnergy >= retentionThreshold, 1, 'first');

    center_oar{it} = center_tmp;
    U_oar{it} = sparse(U(:, 1:k));
    S_oar{it} = sparse(S(1:k, 1:k));
    V_oar{it} = sparse(V(:, 1:k));
end

for it = 1:numel(OARSubIx)
    dij_interval.center(OARSubIx(it), :) = center_oar{it};
    dij_interval.U{it} = U_oar{it};
    dij_interval.S{it} = S_oar{it};
    dij_interval.V{it} = V_oar{it};
end

end
