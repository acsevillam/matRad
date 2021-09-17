function matRad_plotDoseHistos(dij,pln,cst,ixStructure,w)
%MATRAD_CALMAX Summary of this function goes here
%   Detailed explanation goes here

matRad_cfg =  MatRad_Config.instance();

% initialize waitbar
figureWait = waitbar(0,'calculate dose interval for each voxel and bixel...');
% show busy state
set(figureWait,'pointer','watch');

%Set up parent export folder and full file path
if ~(isfolder('dijHistoExport'))
    mkdir(matRad_cfg.matRadRoot, 'dijHistoExport');
end

edges = [0:0.25:10];
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
    
    folderPath=[matRad_cfg.matRadRoot filesep 'dijHistoExport' filesep];
    
    if nnz(dij_tmp)>0
        f=figure;
        set(gcf, 'renderer', 'painters');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [2.5 2.5]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 2.5 2.5]);
        h = histogram(dij_tmp*w,edges)
        xlabel('Dose [Gy]')
        ylabel('Frecuency [Counts]')
        title(['voxel ID =' num2str(i(it))])
        set(gca,...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',8,...
            'FontName','Times')
        saveas(h,[folderPath num2str(i(it)) '_dose.png']);
        close(f);
    end
    
    clear 'dij_tmp';
    
end

%Close Waitbar
if ishandle(figureWait)
    delete(figureWait);
end

end

