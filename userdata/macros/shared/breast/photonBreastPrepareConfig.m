function prepareConfig = photonBreastPrepareConfig(caseID,planTemplate)
% photonBreastPrepareConfig Build common photon breast prepare settings.

prepareConfig = struct();
prepareConfig.caseID = caseID;
prepareConfig.AcquisitionType = 'dicom';
prepareConfig.hlutFileName = 'matRad_default.hlut';
prepareConfig.description = 'breast';
prepareConfig.plan_template = planTemplate;
prepareConfig.radiationMode = 'photons';
prepareConfig.machine = 'Generic';
prepareConfig.bioModel = 'none';
prepareConfig.plan_beams = '7F';
prepareConfig.dicomMetadata = struct();
prepareConfig.resolution = [3 3 3];
prepareConfig.n_cores = 8;

end
