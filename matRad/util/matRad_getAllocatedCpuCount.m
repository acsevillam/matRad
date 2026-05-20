function [cpuCount, source] = matRad_getAllocatedCpuCount(varargin)
% matRad_getAllocatedCpuCount returns scheduler-allocated CPU count
%
% call
%   cpuCount = matRad_getAllocatedCpuCount()
%
% input (optional Name-Value pairs)
%   varargin:    optional Name-Value pairs
%   environment: struct used instead of process environment variables
%
% output
%   cpuCount:    allocated CPU count, or [] when no scheduler limit is known
%   source:      environment variable used to derive cpuCount
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addParameter('environment', [], @(x) isempty(x) || isstruct(x));
p.parse(varargin{:});

environment = p.Results.environment;
cpuCount = [];
source = '';

[cpuCount, source] = matRad_firstPositiveIntegerEnvironmentValue( ...
    environment, {'SLURM_CPUS_PER_TASK', 'SLURM_CPUS_ON_NODE', 'SLURM_NCPUS'});
if ~isempty(cpuCount)
    return
end

jobCpusPerNode = matRad_getEnvironmentValue(environment, 'SLURM_JOB_CPUS_PER_NODE');
[cpuCount, source] = matRad_parseSlurmJobCpusPerNode(jobCpusPerNode);
end

function [value, source] = matRad_firstPositiveIntegerEnvironmentValue(environment, names)
value = [];
source = '';
for i = 1:numel(names)
    candidate = str2double(strtrim(matRad_getEnvironmentValue(environment, names{i})));
    if isfinite(candidate) && candidate >= 1
        value = floor(candidate);
        source = names{i};
        return
    end
end
end

function [cpuCount, source] = matRad_parseSlurmJobCpusPerNode(value)
cpuCount = [];
source = '';
if isempty(value)
    return
end

value = strtrim(char(string(value)));
match = regexp(value, '^([0-9]+)', 'tokens', 'once');
if isempty(match)
    return
end

cpuCount = str2double(match{1});
if ~isfinite(cpuCount) || cpuCount < 1
    cpuCount = [];
    return
end
cpuCount = floor(cpuCount);
source = 'SLURM_JOB_CPUS_PER_NODE';
end

function value = matRad_getEnvironmentValue(environment, name)
if ~isempty(environment)
    if isfield(environment, name)
        value = environment.(name);
    else
        value = '';
    end
else
    value = getenv(name);
end

if isstring(value) && isscalar(value)
    value = char(value);
elseif isnumeric(value) && isscalar(value)
    value = num2str(value);
elseif ~ischar(value)
    value = '';
end
end
