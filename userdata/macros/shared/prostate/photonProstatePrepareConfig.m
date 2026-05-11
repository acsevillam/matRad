function prepareConfig = photonProstatePrepareConfig(prepareConfig)
% photonProstatePrepareConfig Apply photon prostate modality settings.

prepareConfig = ensureProstatePrepareConfig(prepareConfig);
prepareConfig.radiationMode = 'photons';
prepareConfig.machine = 'Generic';
prepareConfig.bioModel = 'none';
prepareConfig.quantityOpt = 'physicalDose';
prepareConfig.plan_beams = '9F';

end
