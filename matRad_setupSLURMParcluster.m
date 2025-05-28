% === Safe parallel pool initialization with custom SLURM-compatible profile ===

function cluster = matRad_setupSLURMParcluster()

    % Define the expected location of the slurm_local.settings profile file
    profilePath = fullfile('matlab_profiles', 'slurm_local.settings');

    % Import the profile only if it is not already available
    if ~any(strcmp(parallel.clusterProfiles, 'slurm_local'))
        if exist(profilePath, 'file')
            parallel.importProfile(profilePath);
        else
            warning('[matRad_sampling] SLURM profile file not found at %s. Using fallback.', profilePath);
        end
    end

    % Create a cluster object from the imported profile
    cluster = parcluster('slurm_local');

    % Configure the job storage location (relevant for SLURM temporary directories)
    if isfolder(getenv('TMPDIR'))
        cluster.JobStorageLocation = getenv('TMPDIR');
    else
        cluster.JobStorageLocation = fullfile(getenv('HOME'), 'matlab_jobs');
        if ~isfolder(cluster.JobStorageLocation)
            mkdir(cluster.JobStorageLocation);
        end
    end

    % Get number of workers from SLURM environment or fallback for local execution
    nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
    if isnan(nWorkers) || nWorkers <= 0
        nWorkers = max(1, feature('numcores') - 2);
    end
    cluster.NumWorkers = nWorkers;

    % Start parallel pool only if it is not already running
    if isempty(gcp('nocreate'))
        fprintf('[matRad_sampling] Starting parpool with %d workers using profile ''slurm_local''...\n', nWorkers);
        parpool(cluster, nWorkers);
    else
        fprintf('[matRad_sampling] parpool already active with %d workers.\n', gcp().NumWorkers);
    end


end
