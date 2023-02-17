function [dij_dummy, pln_dummy,dij,pln,dij_interval] = matRad_calcDoseInterval3(ct,cst,stf,pln,dij,targetStructSel,OARStructSel,k)

matRad_cfg =  MatRad_Config.instance();
[env,envver] = matRad_getEnvironment();

pln_dummy=pln;

% retrieve 2 dummy case scenarios for dose calculation and optimziation
pln_dummy.multScen = matRad_multScen(ct,'nomScen');

% calculate dummy case dij to save interval
dij_dummy = matRad_calcPhotonDose(ct,stf,pln_dummy,cst);

dij_interval.center=sparse(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
dij_interval.radius=sparse(dij.totalNumOfBixels,dij.totalNumOfBixels);

% initialize waitbar
figureWait = waitbar(0,'calculate dose interval for each voxel and bixel...');
% show busy state
set(figureWait,'pointer','watch');

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

for it=1:numel(targetSubIx)

    % Display progress and update text only 200 times
    if matRad_cfg.logLevel > 1
        % Display progress and update text only 200 times
        if mod(it,max(1,round(numel(targetSubIx)/200))) == 0
            matRad_progress(it/max(1,round(numel(targetSubIx)/200)),...
                floor(numel(targetSubIx)/max(1,round(numel(targetSubIx)/200))));
        end

        % update waitbar only 100 times if it is not closed
        if mod(it,round(numel(targetSubIx)/100)) == 0 && ishandle(figureWait)
            waitbar(it/numel(targetSubIx),figureWait);
        end
    end

    scenIx = find(pln.multScen.scenMask);
    dij_tmp=cell2mat(cellfun(@(data) data(targetSubIx(it),:),dij.physicalDose(scenIx),'UniformOutput',false));

    % Interval center dose influence matrix
    dij_interval.center(targetSubIx(it),:)=sum(dij_tmp'*diag(pln.multScen.scenProb),2); % mean(dij_tmp,1);

    % Interval radius dose influence matrix
    dij_interval.radius=dij_interval.radius+(dij_tmp'*diag(pln.multScen.scenProb)*dij_tmp);

    clear 'dij_tmp' 'dij_interval_tmp';

end

dij_interval.U=cell(size(OARSubIx));
dij_interval.S=cell(size(OARSubIx));
dij_interval.V=cell(size(OARSubIx));
dij_interval.OARSubIx=OARSubIx;

for it=1:numel(OARSubIx)

    % Display progress and update text only 200 times
    if matRad_cfg.logLevel > 1
        % Display progress and update text only 200 times
        if mod(it,max(1,round(numel(OARSubIx)/200))) == 0
            matRad_progress(it/max(1,round(numel(OARSubIx)/200)),...
                floor(numel(OARSubIx)/max(1,round(numel(OARSubIx)/200))));
        end

        % update waitbar only 100 times if it is not closed
        if mod(it,round(numel(OARSubIx)/100)) == 0 && ishandle(figureWait)
            waitbar(it/numel(OARSubIx),figureWait);
        end
    end

    scenIx = find(pln.multScen.scenMask);
    dij_tmp=cell2mat(cellfun(@(data) data(OARSubIx(it),:),dij.physicalDose(scenIx),'UniformOutput',false));

    % Interval center dose influence matrix
    dij_interval.center(OARSubIx(it),:)=sum(dij_tmp'*diag(pln.multScen.scenProb),2); % mean(dij_tmp,1);

    % Interval radius dose influence matrix
    dij_interval_tmp.radius=(dij_tmp'*diag(pln.multScen.scenProb)*dij_tmp-dij_interval.center(OARSubIx(it),:)'*dij_interval.center(OARSubIx(it),:));
    [dij_interval.U{it},dij_interval.S{it},dij_interval.V{it}]= svds(dij_interval_tmp.radius,k,'largest');

    clear 'dij_tmp' 'dij_interval_tmp';

end

%Close Waitbar
if ishandle(figureWait)
    delete(figureWait);
end

end