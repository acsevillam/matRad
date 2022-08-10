function [dij_dummy, pln_dummy,dij_interval] = matRad_calcDoseInterval2(ct,cst,stf,pln,dij,structSel)

matRad_cfg =  MatRad_Config.instance();
[env,envver] = matRad_getEnvironment();

pln_dummy=pln;

% retrieve 2 dummy case scenarios for dose calculation and optimziation
pln_dummy.multScen = matRad_multScen(ct,'nomScen');

% calculate dummy case dij to save interval
dij_dummy = matRad_calcPhotonDose(ct,stf,pln_dummy,cst);

dij_interval.center=sparse(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
dij_interval.radius=sparse(dij.totalNumOfBixels,dij.totalNumOfBixels);
%ndSparse(dij_interval.doseGrid.numOfVoxels,dij_interval.totalNumOfBixels,dij_interval.totalNumOfBixels);

% initialize waitbar
figureWait = waitbar(0,'calculate dose interval for each voxel and bixel...');
% show busy state
set(figureWait,'pointer','watch');

% resizing cst to dose cube resolution
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,  dij.ctGrid.y,  dij.ctGrid.z,...
    dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

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

for it=1:numel(subIx)
    
    % Display progress and update text only 200 times
    if matRad_cfg.logLevel > 1
        % Display progress and update text only 200 times
        if mod(it,max(1,round(numel(subIx)/200))) == 0
            matRad_progress(it/max(1,round(numel(subIx)/200)),...
                floor(numel(subIx)/max(1,round(numel(subIx)/200))));
        end
        
        % update waitbar only 100 times if it is not closed
        if mod(it,round(numel(subIx)/100)) == 0 && ishandle(figureWait)
            waitbar(it/numel(subIx),figureWait);
        end
    end
    
    scenIx = find(pln.multScen.scenMask);
    
    dij_tmp=cell2mat(cellfun(@(data) data(subIx(it),:),dij.physicalDose(scenIx),'UniformOutput',false));
    
    % Interval center dose influence matrix 
    dij_interval.center(subIx(it),:)=sum(dij_tmp'*diag(pln.multScen.scenProb),2); % mean(dij_tmp,1);
    
    % Interval radius dose influence matrix
    dij_interval.radius=dij_interval.radius+sparse(dij_tmp'*diag(pln.multScen.scenProb)*dij_tmp);
    
    clear 'dij_tmp';
    
end

%Close Waitbar
if ishandle(figureWait)
    delete(figureWait);
end

end