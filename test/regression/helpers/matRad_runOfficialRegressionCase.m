function result = matRad_runOfficialRegressionCase(targetRoot,caseName)
% Run a regression case in a clean MATLAB/Octave process for targetRoot.

[runner,runnerFile] = writeRunner(targetRoot,caseName);
cleanupObj = onCleanup(@() cleanupFiles(runnerFile,runner.outputFile));

[status,cmdout] = system(runner.command);
if status ~= 0
    error('matRad:RegressionCaseFailed', ...
        'Regression case "%s" failed for root "%s":\n%s', ...
        caseName,targetRoot,cmdout);
end

if exist(runner.outputFile,'file') ~= 2
    error('matRad:RegressionCaseFailed', ...
        'Regression case "%s" did not create output file "%s".', ...
        caseName,runner.outputFile);
end

loaded = load(runner.outputFile,'result');
result = loaded.result;

clear cleanupObj;
cleanupFiles(runnerFile,runner.outputFile);

end

function [runner,runnerFile] = writeRunner(targetRoot,caseName)
    helperDir = fileparts(mfilename('fullpath'));
    runner.outputFile = [tempname() '.mat'];
    runnerFile = [tempname() '.m'];

    fid = fopen(runnerFile,'w');
    if fid < 0
        error('matRad:RegressionCaseFailed', ...
            'Could not create temporary runner "%s".',runnerFile);
    end
    try
        fprintf(fid,'addpath(''%s'');\n',matlabLiteral(helperDir));
        fprintf(fid,'matRad_officialRegressionWorker(''%s'',''%s'',''%s'');\n', ...
            matlabLiteral(targetRoot),matlabLiteral(caseName), ...
            matlabLiteral(runner.outputFile));
        fclose(fid);
    catch ME
        fclose(fid);
        rethrow(ME);
    end

    runner.command = buildExecutionCommand(runnerFile);
end

function command = buildExecutionCommand(runnerFile)
    if exist('OCTAVE_VERSION','builtin')
        executable = getenv('OCTAVE_EXECUTABLE');
        if isempty(executable)
            executable = getenv('OCTAVE');
        end
        if isempty(executable)
            executable = 'octave';
        end
        command = sprintf('"%s" --quiet "%s"',executable,runnerFile);
    else
        executable = fullfile(matlabroot,'bin','matlab');
        if ispc
            executable = [executable '.exe'];
        end
        if exist(executable,'file') ~= 2
            moxunit_throw_test_skipped_exception(sprintf( ...
                'MATLAB executable not found at %s.',executable));
        end
        command = sprintf('"%s" -batch "run(''%s'')"',executable,runnerFile);
    end
end

function value = matlabLiteral(value)
    value = strrep(char(value),'''','''''');
end

function cleanupFiles(varargin)
    for i = 1:numel(varargin)
        if exist(varargin{i},'file') == 2
            delete(varargin{i});
        end
    end
end
