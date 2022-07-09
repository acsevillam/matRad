function [ct,cst] = matRad_loadGeometry(run_config)

    patient=run_config.description;
    caseID=run_config.caseID;
    load([caseID '.mat']);

end