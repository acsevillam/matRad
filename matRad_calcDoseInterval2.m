function [dij_dummy, pln_dummy,dij,pln,dij_interval] = matRad_calcDoseInterval2(ct,cst,stf,pln,dij,structSel)

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
cst_resized = matRad_resizeCstToGrid(cst,dij.ctGrid.x,  dij.ctGrid.y,  dij.ctGrid.z,...
    dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

V = [];
% define voxels for sampling
if ~exist('structSel', 'var') || sum(size(structSel)) == 0
    V = [cst_resized{:,4}];
    % final voxel subset for sampling
    subIx = unique(vertcat(V{:}));
else
    for i=1:size(cst_resized,1)
        for j = 1:numel(structSel)
            if strcmp(structSel{j}, cst_resized{i,2})
                V = [V cst_resized{i,4}{1}];
            end
        end
    end
    % final voxel subset for sampling
    subIx = V;
end

% Calculate scenarios probabilities
scenIx = find(pln.multScen.scenMask);
vProb = zeros(numel(scenIx),1);

for l = 1:numel(scenIx)

    [ctScen,shiftScen,RangeScen] = deal(pln.multScen.linearMask(l,1),pln.multScen.linearMask(l,2),pln.multScen.linearMask(l,3));
    shiftScenMask = find(squeeze(pln.multScen.scenMask(1,:,:)));
    indProb = sub2ind([pln.multScen.totNumShiftScen pln.multScen.totNumRangeScen],shiftScen,RangeScen);

    numCtScen = nnz(pln.multScen.scenMask(:,shiftScen,RangeScen));
    if(numCtScen>1)
        phaseProb=ones(1,numCtScen)/numCtScen;
        vProb(l)=pln.multScen.scenProb(find(shiftScenMask==indProb))/phaseProb(ctScen);
    else
        vProb(l)=pln.multScen.scenProb(find(shiftScenMask==indProb));
    end

end

        if isfield(ct,'dvf')

            if ~isequal(dij.doseGrid.dimensions,ct.dvfDim)  || ~isequal(ct.dvfType,'push')
                matRad_cfg.dispWarning('Dose cube and deformation vector field dimensions are not equal. \n');

                % Instantiate elastic registration
                metadata.nItera = 100;
                metadata.dvfType = 'push';
                register = matRad_ElasticImageRegistration(ct,cst,1,metadata);
                clear 'metadata';

                % Resize dvf according to dose grid
                [ct] = register.calcDVF_resized(dij.doseGrid.resolution);

            end

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

    subIx_deformed = cell(numel(scenIx),1);

    %figure;

    for l = 1:numel(scenIx)

        [ctScen,~,~] = deal(pln.multScen.linearMask(l,1),pln.multScen.linearMask(l,2),pln.multScen.linearMask(l,3));
        [x,y,z] = ind2sub(dij.doseGrid.dimensions,subIx(it));

        tmpCube=zeros(dij.doseGrid.dimensions);
        tmpCube(x,y,z)=1;

        if isfield(ct,'dvf')
            tmpCube_deformed=imwarp(tmpCube, permute(ct.dvf{ctScen},[2 3 4 1]));
            subIx_deformed{l}=find(imwarp(tmpCube, permute(ct.dvf{ctScen},[2 3 4 1])));
        else
            subIx_deformed{l}=subIx;
        end
    
        %{
        nx = ct.cubeDim(1);
        ny = ct.cubeDim(2);
        nz = ct.cubeDim(3);
        [X, Y, Z] = meshgrid(1:size(tmpCube,1),1:size(tmpCube,2), 1:size(tmpCube,3));
        [xq, yq, zq] = meshgrid(linspace(1, size(tmpCube, 2), ny),linspace(1, size(tmpCube, 1), nx), linspace(1, size(tmpCube, 3), nz));
        tmpCube_deformed = interp3(X, Y, Z, tmpCube_deformed, xq, yq, zq,'nearest');

        zDim = ct.dvfDim(3);
        zCtDim = ct.cubeDim(3);
        slice=round(z*zCtDim/zDim);

        subplot(1,numel(scenIx),ctScen);
        matRad_plotSliceWrapper(gca,ct,cst,l,tmpCube_deformed,3,slice,[],[],[],[],[0 1.01]);
        %}
    end



    %dij_tmp=cell2mat(cellfun(@(data) data(subIx(it),:),dij.physicalDose(scenIx),'UniformOutput',false));
    dij_tmp=cell2mat(cellfun(@(data,index) mean(data(index,:),1),dij.physicalDose(scenIx),subIx_deformed,'UniformOutput',false));

    % Interval center dose influence matrix
    dij_interval.center(subIx(it),:)=sum(dij_tmp'*diag(vProb),2); % mean(dij_tmp,1);

    % Interval radius dose influence matrix
    dij_interval.radius=dij_interval.radius+sparse(dij_tmp'*diag(vProb)*dij_tmp);

    clear 'dij_tmp';

end

%Close Waitbar
if ishandle(figureWait)
    delete(figureWait);
end

end