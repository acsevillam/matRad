function test_suite = test_optimizerIPOPT 

test_functions=localfunctions();

initTestSuite;

function test_optimizer_ipopt_construct

    opti = matRad_OptimizerIPOPT();
    assertTrue(isobject(opti));
    assertTrue(isa(opti, 'matRad_OptimizerIPOPT'));
    
function test_optimizer_ipopt_available

    opti = matRad_OptimizerIPOPT();
    assertTrue(opti.IsAvailable());
    assertTrue(matRad_OptimizerIPOPT.IsAvailable()); %Check static
    
function test_optimizer_ipopt_getStatus

    opti = matRad_OptimizerIPOPT();
    [statusmsg, statusflag] = opti.GetStatus();
    assertEqual(statusmsg, 'No Last IPOPT Status Available!');
    assertEqual(statusflag, -1);
        
    % TODO: test other status

% TODO: test optimize function

function test_optimizer_ipopt_apply_options

    opti = matRad_OptimizerIPOPT();
    opti = opti.applyOptions(struct('max_iter',3,'max_cpu_time',12));
    assertEqual(opti.options.max_iter,3);
    assertEqual(opti.options.max_cpu_time,12);

function test_optimizer_ipopt_rejects_unknown_option

    opti = matRad_OptimizerIPOPT();
    assertExceptionThrown(@() opti.applyOptions(struct('definitelyUnknownOption',1)), ...
        'matRad:Error');

function test_optimizer_ipopt_disables_plot_when_gui_disabled

    matRad_cfg = MatRad_Config.instance();
    previousDisableGUI = matRad_cfg.disableGUI;
    cleanup = onCleanup(@() restoreDisableGui(matRad_cfg,previousDisableGUI));
    matRad_cfg.disableGUI = true;

    opti = matRad_OptimizerIPOPT();

    assertFalse(opti.showPlot);
    assertFalse(matRad_isInteractiveSession('requireFigureWindows',true));
    clear cleanup;

function test_create_optimizer_rejects_unknown_optimizer

    assertExceptionThrown(@() matRad_createOptimizer('definitelyUnknownOptimizer'), ...
        'matRad:Error');

function test_create_optimizer_falls_back_to_ipopt_for_unknown_optimizer

    matRad_cfg = MatRad_Config.instance();
    previousKeepLog = matRad_cfg.keepLog;
    previousMessageLog = matRad_cfg.messageLog;
    cleanup = onCleanup(@() restoreKeepLog(matRad_cfg,previousKeepLog, ...
        previousMessageLog));
    matRad_cfg.keepLog = true;
    startLogCount = size(previousMessageLog,1);

    opti = matRad_createOptimizer('definitelyUnknownOptimizer', ...
        'fallbackOptimizer','IPOPT');

    assertTrue(isa(opti,'matRad_OptimizerIPOPT'));
    messages = matRad_cfg.messageLog((startLogCount+1):end,2);
    assertTrue(anyMessageContains(messages,'Fallback to IPOPT'));
    clear cleanup;

function test_create_optimizer_fallback_result_accepts_options

    opti = matRad_createOptimizer('definitelyUnknownOptimizer', ...
        'fallbackOptimizer','IPOPT');
    opti = opti.applyOptions(struct('max_iter',4));

    assertEqual(opti.options.max_iter,4);

function restoreDisableGui(matRad_cfg,previousDisableGUI)
    matRad_cfg.disableGUI = previousDisableGUI;

function restoreKeepLog(matRad_cfg,previousKeepLog,previousMessageLog)
    matRad_cfg.keepLog = previousKeepLog;
    matRad_cfg.messageLog = previousMessageLog;

function tf = anyMessageContains(messages,needle)
    tf = any(cellfun(@(msg) ~isempty(strfind(msg,needle)),messages));
