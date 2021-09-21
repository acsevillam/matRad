function [dij_interval, pln_interval] = matRad_calcDoseInterval2(dij,pln,cst,ixStructure,minPrctile,maxPrctile)
%MATRAD_CALMAX Summary of this function goes here
%   Detailed explanation goes here

dij_interval=dij;
dij_interval.physicalDose=[];
dij_interval.physicalDose{1}=sparse(dij_interval.doseGrid.numOfVoxels,dij_interval.totalNumOfBixels);
dij_interval.physicalDose{2}=sparse(dij_interval.doseGrid.numOfVoxels,dij_interval.totalNumOfBixels);

pln_interval=pln;

matRad_cfg =  MatRad_Config.instance();

% initialize waitbar
figureWait = waitbar(0,'calculate dose interval for each voxel and bixel...');
% show busy state
set(figureWait,'pointer','watch');
i=cst{ixStructure,4}{1,1};

for it = 1:length(i)
    
    % Display progress and update text only 200 times
    if matRad_cfg.logLevel > 1
        % Display progress and update text only 200 times
        if mod(it,max(1,round(length(i)/200))) == 0
            matRad_progress(it/max(1,round(length(i)/200)),...
                floor(length(i)/max(1,round(length(i)/200))));
        end
        
        % update waitbar only 100 times if it is not closed
        if mod(it,round(length(i)/100)) == 0 && ishandle(figureWait)
            waitbar(it/length(i),figureWait);
        end
    end
    
    for ctScen = 1:pln.multScen.numOfCtScen
        for shiftScen = 1:pln.multScen.totNumShiftScen
            for rangeShiftScen = 1:pln.multScen.totNumRangeScen
                if ~isempty(dij.physicalDose{ctScen,shiftScen,rangeShiftScen})
                    if exist('dij_tmp','var') && ~isempty(dij_tmp)
                        dij_tmp=cat(1,dij_tmp,dij.physicalDose{ctScen,shiftScen,rangeShiftScen}(i(it),:));
                    else
                        dij_tmp=dij.physicalDose{ctScen,shiftScen,rangeShiftScen}(i(it),:);
                    end
                end
            end
        end
    end
        
    %Percentile Dose
    dij_interval.physicalDose{1}(i(it),:)=prctile(dij_tmp,minPrctile,1);
    dij_interval.physicalDose{2}(i(it),:)=prctile(dij_tmp,maxPrctile,1);
    
    clear 'dij_tmp';
end

%Close Waitbar
if ishandle(figureWait)
    delete(figureWait);
end

end

