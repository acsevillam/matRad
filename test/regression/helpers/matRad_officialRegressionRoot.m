function officialRoot = matRad_officialRegressionRoot()
% Locate the sibling checkout used as the official regression baseline.

officialRoot = getenv('MATRAD_OFFICIAL_DEV_VAR_RBE_ROOT');
if isempty(officialRoot)
    localRoot = matRad_localRegressionRoot();
    officialRoot = fullfile(fileparts(localRoot),'dev_varRBErobOpt_official');
end

officialRoot = char(officialRoot);

end
