function ok = matRad_testSLURMParpool()
% matRad_testSLURMParpool checks if the SLURM-compatible parallel pool can be initialized and is operational.
%
% Output:
%   ok - true if everything works correctly, false otherwise

ok = false;

try
    fprintf('[matRad_testSLURMParpool] Starting test...\n');

    % Initialize the cluster and parpool
    matRad_setupSLURMParcluster();

    % Run a small parallel test
    N = 10;
    A = zeros(1, N);

    fprintf('[matRad_testSLURMParpool] Running parallel test using parfor...\n');
    parfor i = 1:N
        A(i) = i^2;
    end

    % Validate result
    if isequal(A, (1:N).^2)
        fprintf('[matRad_testSLURMParpool] Test successful. parpool is working correctly.\n');
        ok = true;
    else
        warning('[matRad_testSLURMParpool] Parallel execution returned incorrect result.');
    end

    % Optional: close the pool after test
    delete(gcp('nocreate'));

catch ME
    warning('[matRad_testSLURMParpool] Test failed: %s', ME.message);
end
end
