function [dij_dummy, pln_dummy,dij,pln,dij_interval] = matRad_calcDoseInterval3b(ct,cst,stf,pln,dij,targetStructSel,OARStructSel,k)

matRad_cfg =  MatRad_Config.instance();
[env,envver] = matRad_getEnvironment();

pln_dummy=pln;

% retrieve 2 dummy case scenarios for dose calculation and optimziation
pln_dummy.multScen = matRad_multScen(ct,'nomScen');

% calculate dummy case dij to save interval
dij_dummy = matRad_calcPhotonDose(ct,stf,pln_dummy,cst);

scenIx = find(pln.multScen.scenMask);
dij_list=dij.physicalDose(scenIx);
scenProb=pln.multScen.scenProb;

% resizing cst to dose cube resolution
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,  dij.ctGrid.y,  dij.ctGrid.z,...
    dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

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

p = gcp(); % If no pool, create new one.

columnHeadings = {'Ix' 'center','radius'};
% Preallocate structure
dij_interval_target(1:numel(targetSubIx)) = cell2struct(repmat({[]},numel(columnHeadings),1),columnHeadings,1);

if exist('parfor_progress', 'file') == 2
    FlagParforProgressDisp = true;
    parfor_progress(round(numel(targetSubIx)/1000));  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
else
    matRad_cfg.dispInfo('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n');
    FlagParforProgressDisp = false;
end

parfor it=1:numel(targetSubIx)
    
    Ix=targetSubIx(it);
    dij_tmp=cell2mat(cellfun(@(data) data(Ix,:),dij_list,'UniformOutput',false));

    % Voxel index
    dij_interval_target(it).Ix=Ix;

    % Interval center dose influence matrix
    dij_interval_target(it).center=sum(dij_tmp'*diag(scenProb),2); % mean(dij_tmp,1);

    % Interval radius dose influence matrix
    dij_interval_target(it).radius=(dij_tmp'*diag(scenProb)*dij_tmp);

    if FlagParforProgressDisp && mod(it,1000)==0
        parfor_progress;
    end
end

if FlagParforProgressDisp
    parfor_progress(0);
end

dij_interval.center=sparse(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
dij_interval.radius=sparse(dij.totalNumOfBixels,dij.totalNumOfBixels);
dij_interval.targetSubIx=targetSubIx;

for it=1:numel(targetSubIx)
    dij_interval.center(targetSubIx(it),:)=dij_interval_target(it).center;
    % Interval radius dose influence matrix
    dij_interval.radius=dij_interval.radius+dij_interval_target(it).radius;
end

clear 'dij_interval_target';

columnHeadings = {'Ix' 'center','U','S','V'};
% Preallocate structure
dij_interval_OAR(1:numel(OARSubIx)) = cell2struct(repmat({[]},numel(columnHeadings),1),columnHeadings,1);

if exist('parfor_progress', 'file') == 2
    FlagParforProgressDisp = true;
    parfor_progress(round(numel(OARSubIx)/1000));  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
else
    matRad_cfg.dispInfo('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n');
    FlagParforProgressDisp = false;
end

parfor it=1:numel(OARSubIx)

    Ix=OARSubIx(it);
    dij_tmp=cell2mat(cellfun(@(data) data(Ix,:),dij_list,'UniformOutput',false));

    % Voxel index
    dij_interval_OAR(it).Ix=Ix;

    % Interval center dose influence matrix
    dij_interval_OAR(it).center=sum(dij_tmp'*diag(scenProb),2); % mean(dij_tmp,1);

    % Interval radius dose influence matrix
    radius_tmp=(dij_tmp'*diag(scenProb)*dij_tmp-dij_interval_OAR(it).center'*dij_interval_OAR(it).center);
    [dij_interval_OAR(it).U,dij_interval_OAR(it).S,dij_interval_OAR(it).V]= svds(radius_tmp,k,'largest'); 

    if FlagParforProgressDisp && mod(it,1000)==0
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