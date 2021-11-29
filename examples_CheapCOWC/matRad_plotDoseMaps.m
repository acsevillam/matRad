
%% Clear variables
clear;
clc;

%% set matRad runtime configuration
matRad_rc
param.logLevel=1;

%% Create an interactive plot to slide through axial slices

quantityOpt='physicalDose';
radiationMode="photons";

if ~exist('rootPath','var') || isempty(rootPath)
    run_config.rootPath = matRad_cfg.matRadRoot;
else
    run_config.rootPath = rootPath;
end

if run_config.robustness == "c-COWC"
    output_folder = ['output' filesep filesep run_config.description filesep run_config.robustness filesep run_config.plan_target filesep run_config.plan_beams filesep run_config.plan_objectives filesep run_config.scen_mode filesep num2str(run_config.wcFactor) filesep num2str(run_config.beta1) '_to_' num2str(run_config.beta2) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
else
    output_folder = ['output' filesep filesep run_config.description filesep run_config.robustness filesep run_config.plan_target filesep run_config.plan_beams filesep run_config.plan_objectives filesep run_config.scen_mode filesep num2str(run_config.wcFactor) filesep datestr(datetime,'yyyy-mm-dd HH-MM-SS')];
end

%Set up parent export folder and full file path
if ~(isfolder(output_folder))
    mkdir(run_config.rootPath, output_folder);
end

folderPath = [run_config.rootPath filesep output_folder];

%% Plot beam doses

run_config.doseWindow(2)=run_config.doseWindow(2)/2;

if radiationMode == "photons" 
    for numBeam=1:pln_robust.propStf.numOfBeams

        if(pln_robust.propStf.numOfBeams>1)
            quantityMap=[quantityOpt '_beam' num2str(round(numBeam))];
        else
            quantityMap=quantityOpt;
        end

        plane      = 3;
        doseWindow = [0 max([run_config.doseWindow(2)])];
        doseIsoLevels = 0; %linspace(0.1 * maxDose,maxDose,10);
        f = figure;
        title([quantityMap]);
        set(gcf,'position',[10,10,550,400]);
        numSlices = ct.cubeDim(3);
        slice = round(pln_robust.propStf.isoCenter(1,3)./ct.resolution.z);
        matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUI_robust.(quantityMap)*pln_robust.numOfFractions,plane,slice,[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]');
        b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
              'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
        b.Callback = @(es,ed)  matRad_plotSliceWrapper(gca,ct,cst_robust,1,resultGUI_robust.(quantityMap)*pln_robust.numOfFractions,plane,round(es.Value),[],[],colorcube,[],doseWindow,doseIsoLevels,[],'Dose [Gy]');

        savefig([folderPath filesep 'dose_robust_' quantityMap '.fig']);   
    end
end