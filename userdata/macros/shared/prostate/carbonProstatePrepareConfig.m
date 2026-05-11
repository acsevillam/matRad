function prepareConfig = carbonProstatePrepareConfig(prepareConfig)
% carbonProstatePrepareConfig Apply carbon prostate modality settings.

prepareConfig = ensureProstatePrepareConfig(prepareConfig);
prepareConfig.radiationMode = 'carbon';
prepareConfig.machine = 'Generic';
prepareConfig.bioModel = 'LEM';
prepareConfig.quantityOpt = 'RBExD';
prepareConfig.plan_beams = '2F';

end
