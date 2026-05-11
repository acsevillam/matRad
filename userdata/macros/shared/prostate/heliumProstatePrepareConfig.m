function prepareConfig = heliumProstatePrepareConfig(prepareConfig)
% heliumProstatePrepareConfig Apply helium prostate modality settings.

prepareConfig = ensureProstatePrepareConfig(prepareConfig);
prepareConfig.radiationMode = 'helium';
prepareConfig.machine = 'Generic';
prepareConfig.bioModel = 'HEL';
prepareConfig.quantityOpt = 'RBExD';
prepareConfig.plan_beams = '2F';

end
