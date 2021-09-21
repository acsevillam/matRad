function [dij_interval, pln_interval] = matRad_calcDoseInterval(dij,pln,structure,minPrctile,maxPrctile)
%MATRAD_CALMAX Summary of this function goes here
%   Detailed explanation goes here

dij_interval=dij;
dij_interval.physicalDose=[];
pln_interval=pln;

matRad_cfg =  MatRad_Config.instance();

% initialize waitbar
figureWait = waitbar(0,'calculate dose interval for each voxel and bixel...');
% show busy state
set(figureWait,'pointer','watch');
structure_size=length(structure);

for j = 1:structure_size
    
    % Display progress and update text only 200 times
    if matRad_cfg.logLevel > 1
        % Display progress and update text only 200 times
        if mod(j,max(1,round(structure_size/200))) == 0
            matRad_progress(j/max(1,round(structure_size/200)),...
                floor(structure_size/max(1,round(structure_size/200))));
        end
        
        % update waitbar only 100 times if it is not closed
        if mod(j,round(structure_size/100)) == 0 && ishandle(figureWait)
            waitbar(j/structure_size,figureWait);
        end
    end
    
    for ctScen = 1:pln.multScen.numOfCtScen
        for shiftScen = 1:pln.multScen.totNumShiftScen
            for rangeShiftScen = 1:pln.multScen.totNumRangeScen
                if ~isempty(dij.physicalDose{ctScen,shiftScen,rangeShiftScen})
                    if exist('dij_tmp','var') && ~isempty(dij_tmp)
                        dij_tmp=cat(2,dij_tmp,dij.physicalDose{ctScen,shiftScen,rangeShiftScen}(:,j));
                    else
                        dij_tmp=dij.physicalDose{ctScen,shiftScen,rangeShiftScen}(:,j);
                    end
                end
            end
        end
    end
    
    if ~isempty(dij_interval.physicalDose)

        %Percentile Dose
        dij_interval.physicalDose{1}=cat(2,dij_interval.physicalDose{1},prctile(dij_tmp,minPrctile,2));
        dij_interval.physicalDose{2}=cat(2,dij_interval.physicalDose{2},prctile(dij_tmp,maxPrctile,2));

    else
        
        %Percentile Dose
        dij_interval.physicalDose{1}=prctile(dij_tmp,minPrctile,2);
        dij_interval.physicalDose{2}=prctile(dij_tmp,maxPrctile,2);
        
    end
    
    clear 'dij_tmp';
end

%Close Waitbar
if ishandle(figureWait)
    delete(figureWait);
end

end

