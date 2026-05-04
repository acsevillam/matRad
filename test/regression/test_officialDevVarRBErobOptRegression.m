function test_suite = test_officialDevVarRBErobOptRegression
% Regression tests against an isolated official dev_varRBErobOpt checkout.

test_functions = localfunctions();

initTestSuite;

function test_scenario_generation_matches_official_dev_varRBErobOpt
    officialRoot = matRad_officialRegressionRoot();
    assertOfficialCheckoutAvailable(officialRoot);

    localResult = matRad_runOfficialRegressionCase(matRad_localRegressionRoot(), ...
        'scenario_generation');
    officialResult = matRad_runOfficialRegressionCase(officialRoot, ...
        'scenario_generation');

    matRad_compareRegressionStructs(officialResult,localResult, ...
        'scenario_generation',1e-10);

function test_dose_engine_initialization_matches_official_dev_varRBErobOpt
    officialRoot = matRad_officialRegressionRoot();
    assertOfficialCheckoutAvailable(officialRoot);

    localResult = matRad_runOfficialRegressionCase(matRad_localRegressionRoot(), ...
        'dose_engine_initialization');
    officialResult = matRad_runOfficialRegressionCase(officialRoot, ...
        'dose_engine_initialization');

    matRad_compareRegressionStructs(officialResult,localResult, ...
        'dose_engine_initialization',1e-10);

function test_pencil_beam_dose_matches_official_dev_varRBErobOpt
    officialRoot = matRad_officialRegressionRoot();
    assertOfficialCheckoutAvailable(officialRoot);

    localResult = matRad_runOfficialRegressionCase(matRad_localRegressionRoot(), ...
        'pencil_beam_dose');
    officialResult = matRad_runOfficialRegressionCase(officialRoot, ...
        'pencil_beam_dose');

    matRad_compareRegressionStructs(officialResult,localResult, ...
        'pencil_beam_dose',1e-8);

function test_particle_pencil_beam_biological_models_match_official_dev_varRBErobOpt
    officialRoot = matRad_officialRegressionRoot();
    assertOfficialCheckoutAvailable(officialRoot);

    localResult = matRad_runOfficialRegressionCase(matRad_localRegressionRoot(), ...
        'particle_bio_pencil_beam_dose');
    officialResult = matRad_runOfficialRegressionCase(officialRoot, ...
        'particle_bio_pencil_beam_dose');

    matRad_compareRegressionStructs(officialResult,localResult, ...
        'particle_bio_pencil_beam_dose',1e-8);

function assertOfficialCheckoutAvailable(officialRoot)
    if exist(officialRoot,'dir') ~= 7 || ...
            exist(fullfile(officialRoot,'matRad_rc.m'),'file') ~= 2
        moxunit_throw_test_skipped_exception(sprintf([ ...
            'Official dev_varRBErobOpt checkout not found at %s. ', ...
            'Set MATRAD_OFFICIAL_DEV_VAR_RBE_ROOT to override.'],officialRoot));
    end
