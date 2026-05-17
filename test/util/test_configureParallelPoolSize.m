function test_suite = test_configureParallelPoolSize

test_functions = localfunctions();

initTestSuite;

function test_configureParallelPoolCreatesAndReusesRequestedSize
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception(['Parallel Computing ', ...
                                          'Toolbox is unavailable.']);
end

cleanup = helper_preserveParallelPool(); %#ok<NASGU>
helper_deleteActiveParallelPool();

pPool = matRad_configureParallelPoolSize(1, 'test pool creation', ...
                                         MatRad_Config.instance());
assertEqual(pPool.NumWorkers, 1);

pPool = matRad_configureParallelPoolSize(1, 'test pool reuse', ...
                                         MatRad_Config.instance());
assertEqual(pPool.NumWorkers, 1);
activePool = gcp('nocreate');
assertEqual(activePool.NumWorkers, 1);

function test_configureParallelPoolIncreasesRequestedSize
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception(['Parallel Computing ', ...
                                          'Toolbox is unavailable.']);
end
if helper_maxParallelPoolWorkers() < 2
    moxunit_throw_test_skipped_exception(['Local parallel pool ', ...
                                          'supports fewer than two workers.']);
end

cleanup = helper_preserveParallelPool(); %#ok<NASGU>

pPool = matRad_configureParallelPoolSize(1, ...
                                         'test pool increase setup', MatRad_Config.instance());
assertEqual(pPool.NumWorkers, 1);

pPool = matRad_configureParallelPoolSize(2, 'test pool increase', ...
                                         MatRad_Config.instance());
assertEqual(pPool.NumWorkers, 2);
activePool = gcp('nocreate');
assertEqual(activePool.NumWorkers, 2);

function test_configureParallelPoolReducesRequestedSize
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception(['Parallel Computing ', ...
                                          'Toolbox is unavailable.']);
end
if helper_maxParallelPoolWorkers() < 2
    moxunit_throw_test_skipped_exception(['Local parallel pool ', ...
                                          'supports fewer than two workers.']);
end

cleanup = helper_preserveParallelPool(); %#ok<NASGU>

pPool = matRad_configureParallelPoolSize(2, ...
                                         'test pool reduction setup', MatRad_Config.instance());
assertEqual(pPool.NumWorkers, 2);

pPool = matRad_configureParallelPoolSize(1, 'test pool reduction', ...
                                         MatRad_Config.instance());
assertEqual(pPool.NumWorkers, 1);
activePool = gcp('nocreate');
assertEqual(activePool.NumWorkers, 1);

function test_configureParallelPoolRejectsInvalidTargetWithoutPoolChange
if helper_parallelComputingAvailable()
    cleanup = helper_preserveParallelPool(); %#ok<NASGU>
    pPool = matRad_configureParallelPoolSize(1, ...
                                             'test invalid target setup', MatRad_Config.instance());
    originalNumWorkers = pPool.NumWorkers;
else
    originalNumWorkers = helper_currentParallelPoolWorkers();
end
matRadCfg = MatRad_Config.instance();

assertExceptionThrown(@() matRad_configureParallelPoolSize(1.5, ...
                                                           'test invalid target', matRadCfg), 'matRad:Error');
assertEqual(helper_currentParallelPoolWorkers(), originalNumWorkers);

assertExceptionThrown(@() matRad_configureParallelPoolSize(0, ...
                                                           'test invalid target', matRadCfg), 'matRad:Error');
assertEqual(helper_currentParallelPoolWorkers(), originalNumWorkers);

function test_configureParallelPoolRestoresOriginalPoolAfterResizeError
if ~helper_parallelComputingAvailable()
    moxunit_throw_test_skipped_exception(['Parallel Computing ', ...
                                          'Toolbox is unavailable.']);
end
maxWorkers = helper_maxParallelPoolWorkers();
if isempty(maxWorkers) || ~isnumeric(maxWorkers) || maxWorkers < 1
    moxunit_throw_test_skipped_exception(['Could not determine a ', ...
                                          'reliable invalid worker count for this environment.']);
end

cleanup = helper_preserveParallelPool(); %#ok<NASGU>

pPool = matRad_configureParallelPoolSize(1, ...
                                         'test pool restore setup', MatRad_Config.instance());
assertEqual(pPool.NumWorkers, 1);

invalidWorkerCount = maxWorkers + 1;
assertExceptionThrown(@() matRad_configureParallelPoolSize( ...
                                                           invalidWorkerCount, 'test pool restore', MatRad_Config.instance()), ...
                      'matRad:Error');

activePool = gcp('nocreate');
assertTrue(~isempty(activePool));
assertEqual(activePool.NumWorkers, 1);

function cleanup = helper_preserveParallelPool()
try
    pPool = gcp('nocreate');
catch
    pPool = [];
end
if isempty(pPool)
    originalNumWorkers = [];
else
    originalNumWorkers = pPool.NumWorkers;
end
cleanup = onCleanup(@() helper_restoreParallelPool(originalNumWorkers));

function helper_restoreParallelPool(originalNumWorkers)
helper_deleteActiveParallelPool();
if ~isempty(originalNumWorkers)
    try
        parpool(originalNumWorkers);
    catch
    end
end

function helper_deleteActiveParallelPool()
try
    pPool = gcp('nocreate');
catch
    return
end
if ~isempty(pPool)
    delete(pPool);
end

function numWorkers = helper_currentParallelPoolWorkers()
numWorkers = [];
if exist('gcp', 'file') ~= 2
    return
end
try
    pPool = gcp('nocreate');
catch
    return
end
if ~isempty(pPool)
    numWorkers = pPool.NumWorkers;
end

function maxWorkers = helper_maxParallelPoolWorkers()
maxWorkers = 0;
try
    cluster = parcluster();
    maxWorkers = cluster.NumWorkers;
catch
end

if isempty(maxWorkers) || ~isnumeric(maxWorkers) || maxWorkers < 1
    try
        maxWorkers = feature('numcores');
    catch
        maxWorkers = 0;
    end
end

function available = helper_parallelComputingAvailable()
available = false;
if exist('parpool', 'file') ~= 2 || exist('gcp', 'file') ~= 2
    return
end

try
    [available, ~] = license('checkout', 'Distrib_Computing_Toolbox');
catch
    available = false;
end

if isempty(available)
    available = false;
end
available = logical(available);
