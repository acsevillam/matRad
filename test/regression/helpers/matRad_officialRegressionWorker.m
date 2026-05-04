function matRad_officialRegressionWorker(targetRoot,caseName,outFile)
% Worker executed in a fresh process for one checkout and one case.

targetRoot = char(targetRoot);
caseName = char(caseName);
outFile = char(outFile);

helperDir = fileparts(mfilename('fullpath'));
regressionDir = fileparts(helperDir);

oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));

if exist(fullfile(targetRoot,'matRad_rc.m'),'file') ~= 2
    error('matRad:RegressionCaseFailed', ...
        'Target root "%s" is not a matRad checkout.',targetRoot);
end

addpath(targetRoot);
cd(targetRoot);
matRad_cfg = matRad_rc(false);
matRad_cfg.setDefaultPropertiesForTesting();
matRad_cfg.disableGUI = true;
matRad_cfg.logLevel = 1;

addpath(regressionDir);
addpath(helperDir);

result = matRad_officialRegressionCase(caseName);

save(outFile,'result');

clear cleanupObj;
cd(oldDir);

end
