function localRoot = matRad_localRegressionRoot()
% Return the matRad root for the checkout running the MOxUnit test.

try
    matRad_cfg = MatRad_Config.instance();
    localRoot = matRad_cfg.matRadRoot;
catch
    helperDir = fileparts(mfilename('fullpath'));
    regressionDir = fileparts(helperDir);
    testDir = fileparts(regressionDir);
    localRoot = fileparts(testDir);
end

localRoot = char(localRoot);

end
