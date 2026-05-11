function prepareConfig = protonProstatePrepareConfig(prepareConfig)
% protonProstatePrepareConfig Apply proton prostate modality settings.

prepareConfig = ensureProstatePrepareConfig(prepareConfig);
prepareConfig.radiationMode = 'protons';
prepareConfig.machine = 'Generic';
prepareConfig.bioModel = 'constRBE';
prepareConfig.quantityOpt = 'RBExD';
prepareConfig.plan_beams = '2F';

end
