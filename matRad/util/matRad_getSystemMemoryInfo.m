function memoryInfo = matRad_getSystemMemoryInfo(varargin)
% matRad_getSystemMemoryInfo returns system memory information in bytes
%
% call
%   memoryInfo = matRad_getSystemMemoryInfo()
%   memoryInfo = matRad_getSystemMemoryInfo('reserveFraction',reserveFraction)
%   memoryInfo = matRad_getSystemMemoryInfo('minReserveBytes',minReserveBytes)
%
% input (optional Name-Value pairs)
%   varargin:        optional Name-Value pairs
%   reserveFraction:  fraction of total system memory kept in reserve
%   minReserveBytes: lower bound for the reserved memory budget
%   includeJobLimits: respect scheduler/cgroup memory limits (default true)
%
% output
%   memoryInfo:       struct with available, total, reserved, and usable
%                     system memory estimates in bytes. The source field
%                     reports which operating-system interface was used.
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
p.addParameter('reserveFraction', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x < 1);
p.addParameter('minReserveBytes', 0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('includeJobLimits', true, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
p.addParameter('environment', [], @(x) isempty(x) || isstruct(x));
p.addParameter('cgroupRoot', '/sys/fs/cgroup', @(x) ischar(x) || (isstring(x) && isscalar(x)));
p.parse(varargin{:});

memoryInfo = struct('availableBytes', [], 'totalBytes', [], ...
                    'reserveBytes', [], 'usableBytes', [], 'source', '');

environment = p.Results.environment;
cgroupRoot = char(p.Results.cgroupRoot);
inSlurm = false;
source = '';
availableBytes = [];
totalBytes = [];

if logical(p.Results.includeJobLimits)
    [availableBytes, totalBytes, source, inSlurm] = ...
        matRad_getJobLimitedMemoryBytes(environment, cgroupRoot);
end

if isempty(availableBytes)
    if inSlurm
        memoryInfo.source = source;
        return
    end
    [availableBytes, totalBytes, source] = matRad_getAvailableSystemMemoryBytes();
end

memoryInfo.availableBytes = availableBytes;
memoryInfo.totalBytes = totalBytes;
memoryInfo.source = source;

if isempty(availableBytes)
    return
end

if isempty(totalBytes)
    reserveBaseBytes = availableBytes;
else
    reserveBaseBytes = totalBytes;
end

memoryInfo.reserveBytes = max(reserveBaseBytes * p.Results.reserveFraction, ...
                              p.Results.minReserveBytes);
memoryInfo.usableBytes = max(0, availableBytes - memoryInfo.reserveBytes);
end

function [availableBytes, totalBytes, source, inSlurm] = matRad_getJobLimitedMemoryBytes(environment, cgroupRoot)
availableBytes = [];
totalBytes = [];
source = '';
inSlurm = matRad_isSlurmJob(environment);

[availableBytes, totalBytes, source] = matRad_getCgroupMemoryBytes(cgroupRoot);
if ~isempty(availableBytes)
    return
end

if inSlurm
    [availableBytes, totalBytes, source] = matRad_getSlurmMemoryBytes(environment);
    if isempty(source)
        source = 'slurm';
    end
end
end

function tf = matRad_isSlurmJob(environment)
tf = false;
slurmFields = {'SLURM_JOB_ID', 'SLURM_JOBID', 'SLURM_STEP_ID', 'SLURM_PROCID'};
for i = 1:numel(slurmFields)
    if ~isempty(matRad_getEnvironmentValue(environment, slurmFields{i}))
        tf = true;
        return
    end
end
end

function [availableBytes, totalBytes, source] = matRad_getSlurmMemoryBytes(environment)
availableBytes = [];
totalBytes = [];
source = '';

memPerNodeBytes = matRad_parseSlurmMemoryBytes( ...
    matRad_getEnvironmentValue(environment, 'SLURM_MEM_PER_NODE'));
if ~isempty(memPerNodeBytes)
    availableBytes = memPerNodeBytes;
    totalBytes = memPerNodeBytes;
    source = 'slurm:SLURM_MEM_PER_NODE';
    return
end

memPerCpuBytes = matRad_parseSlurmMemoryBytes( ...
    matRad_getEnvironmentValue(environment, 'SLURM_MEM_PER_CPU'));
[cpuCount, cpuSource] = matRad_getAllocatedCpuCount('environment', environment);
if ~isempty(memPerCpuBytes) && ~isempty(cpuCount)
    totalBytes = memPerCpuBytes * cpuCount;
    availableBytes = totalBytes;
    source = ['slurm:SLURM_MEM_PER_CPU*' cpuSource];
end
end

function memoryBytes = matRad_parseSlurmMemoryBytes(value)
memoryBytes = [];
if isempty(value)
    return
end

value = strtrim(char(string(value)));
match = regexp(value, '^([0-9]+(?:\.[0-9]+)?)([KkMmGgTt]?)$', 'tokens', 'once');
if isempty(match)
    return
end

number = str2double(match{1});
if ~isfinite(number) || number <= 0
    return
end

switch upper(match{2})
    case 'K'
        memoryBytes = number * 1024;
    case {'', 'M'}
        memoryBytes = number * 1024^2;
    case 'G'
        memoryBytes = number * 1024^3;
    case 'T'
        memoryBytes = number * 1024^4;
end
end

function [availableBytes, totalBytes, source] = matRad_getCgroupMemoryBytes(cgroupRoot)
availableBytes = [];
totalBytes = [];
source = '';

if isempty(cgroupRoot) || exist(cgroupRoot, 'dir') ~= 7
    return
end

candidates = matRad_getCgroupMemoryCandidates(cgroupRoot);
bestTotalBytes = [];
bestAvailableBytes = [];
for i = 1:numel(candidates)
    candidate = candidates{i};
    [limitBytes, currentBytes] = matRad_readCgroupMemoryCandidate(candidate);
    if isempty(limitBytes)
        continue
    end
    if isempty(currentBytes)
        currentBytes = 0;
    end
    candidateAvailableBytes = max(0, limitBytes - currentBytes);
    if isempty(bestTotalBytes) || limitBytes < bestTotalBytes
        bestTotalBytes = limitBytes;
        bestAvailableBytes = candidateAvailableBytes;
    end
end

if ~isempty(bestTotalBytes)
    totalBytes = bestTotalBytes;
    availableBytes = bestAvailableBytes;
    source = 'cgroup';
end
end

function candidates = matRad_getCgroupMemoryCandidates(cgroupRoot)
candidates = {};
candidates{end + 1} = cgroupRoot;
candidates{end + 1} = fullfile(cgroupRoot, 'memory');

try
    cgroupText = fileread('/proc/self/cgroup');
catch
    return
end

lines = regexp(cgroupText, '\r?\n', 'split');
for i = 1:numel(lines)
    line = strtrim(lines{i});
    if isempty(line)
        continue
    end
    fields = regexp(line, ':', 'split');
    if numel(fields) ~= 3
        continue
    end

    cgroupPath = matRad_normalizeCgroupPath(fields{3});
    if strcmp(fields{1}, '0')
        candidates{end + 1} = fullfile(cgroupRoot, cgroupPath); %#ok<AGROW>
    elseif ~isempty(strfind(fields{2}, 'memory'))
        candidates{end + 1} = fullfile(cgroupRoot, 'memory', cgroupPath); %#ok<AGROW>
    end
end
end

function pathValue = matRad_normalizeCgroupPath(pathValue)
pathValue = char(string(pathValue));
while ~isempty(pathValue) && pathValue(1) == '/'
    pathValue = pathValue(2:end);
end
end

function [limitBytes, currentBytes] = matRad_readCgroupMemoryCandidate(candidate)
limitBytes = [];
currentBytes = [];

limitBytes = matRad_readCgroupLimitFile(fullfile(candidate, 'memory.max'));
if isempty(limitBytes)
    limitBytes = matRad_readCgroupLimitFile(fullfile(candidate, 'memory.limit_in_bytes'));
end
if isempty(limitBytes)
    return
end

currentBytes = matRad_readCgroupUsageFile(fullfile(candidate, 'memory.current'));
if isempty(currentBytes)
    currentBytes = matRad_readCgroupUsageFile(fullfile(candidate, 'memory.usage_in_bytes'));
end
end

function valueBytes = matRad_readCgroupLimitFile(filePath)
valueBytes = matRad_readCgroupNumericFile(filePath);
if isempty(valueBytes) || valueBytes <= 0 || valueBytes >= 2^60
    valueBytes = [];
end
end

function valueBytes = matRad_readCgroupUsageFile(filePath)
valueBytes = matRad_readCgroupNumericFile(filePath);
if isempty(valueBytes) || valueBytes < 0
    valueBytes = [];
end
end

function valueBytes = matRad_readCgroupNumericFile(filePath)
valueBytes = [];
try
    valueText = strtrim(fileread(filePath));
catch
    return
end

if strcmp(valueText, 'max')
    return
end
valueBytes = str2double(valueText);
if ~isfinite(valueBytes)
    valueBytes = [];
end
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

function [availableBytes, totalBytes, source] = matRad_getAvailableSystemMemoryBytes()
availableBytes = [];
totalBytes = [];
source = '';

if ispc
    [availableBytes, totalBytes, source] = matRad_getWindowsSystemMemoryBytes();
elseif ismac
    [availableBytes, totalBytes, source] = matRad_getMacSystemMemoryBytes();
elseif isunix
    [availableBytes, totalBytes, source] = matRad_getUnixSystemMemoryBytes();
end
end

function [availableBytes, totalBytes, source] = matRad_getWindowsSystemMemoryBytes()
try
    [~, systemView] = memory();
    availableBytes = double(systemView.PhysicalMemory.Available);
    totalBytes = double(systemView.PhysicalMemory.Total);
    source = 'memory';
    return
catch
end

[availableBytes, totalBytes, source] = matRad_getWindowsWmicMemoryBytes();
end

function [availableBytes, totalBytes, source] = matRad_getWindowsWmicMemoryBytes()
availableBytes = [];
totalBytes = [];
source = '';

[status, wmicText] = system('wmic OS get FreePhysicalMemory,TotalVisibleMemorySize /Value');
if status ~= 0
    return
end

availableKb = matRad_getKeyValueNumber(wmicText, 'FreePhysicalMemory');
totalKb = matRad_getKeyValueNumber(wmicText, 'TotalVisibleMemorySize');

if ~isnan(availableKb)
    availableBytes = availableKb * 1024;
end
if ~isnan(totalKb)
    totalBytes = totalKb * 1024;
end
if ~isempty(availableBytes)
    source = 'wmic';
end
end

function [availableBytes, totalBytes, source] = matRad_getMacSystemMemoryBytes()
availableBytes = [];
totalBytes = [];
source = '';

[memoryPressureAvailableBytes, memoryPressureTotalBytes] = matRad_getMacMemoryPressureBytes();
[vmStatAvailableBytes, vmStatTotalBytes] = matRad_getMacVmStatMemoryBytes();

availableCandidates = [memoryPressureAvailableBytes vmStatAvailableBytes];
availableCandidates = availableCandidates(~isnan(availableCandidates));
if isempty(availableCandidates)
    return
end

availableBytes = min(availableCandidates);

totalCandidates = [memoryPressureTotalBytes vmStatTotalBytes];
totalCandidates = totalCandidates(~isnan(totalCandidates));
if ~isempty(totalCandidates)
    totalBytes = max(totalCandidates);
end

if ~isnan(vmStatAvailableBytes) && availableBytes == vmStatAvailableBytes
    source = 'vm_stat';
else
    source = 'memory_pressure';
end

if ~isnan(memoryPressureAvailableBytes) && ~isnan(vmStatAvailableBytes)
    source = 'memory_pressure/vm_stat';
end
end

function [availableBytes, totalBytes] = matRad_getMacMemoryPressureBytes()
availableBytes = NaN;
totalBytes = NaN;

[status, memoryPressureText] = system('memory_pressure');
if status == 0
    totalMatch = regexp(memoryPressureText, 'The system has\s+(\d+)\s+\(', 'tokens', 'once');
    freePercentMatch = regexp(memoryPressureText, 'System-wide memory free percentage:\s+([0-9.]+)%', 'tokens', 'once');
    if ~isempty(totalMatch) && ~isempty(freePercentMatch)
        totalBytes = str2double(totalMatch{1});
        freePercent = str2double(freePercentMatch{1});
        if isfinite(totalBytes) && isfinite(freePercent)
            availableBytes = totalBytes * freePercent / 100;
            return
        end
    end
end
end

function [availableBytes, totalBytes] = matRad_getMacVmStatMemoryBytes()
availableBytes = NaN;
totalBytes = NaN;

[status, vmStatText] = system('vm_stat');
if status ~= 0
    return
end

pageSizeMatch = regexp(vmStatText, 'page size of\s+(\d+)\s+bytes', 'tokens', 'once');
if isempty(pageSizeMatch)
    pageSize = 4096;
else
    pageSize = str2double(pageSizeMatch{1});
end

freePages = matRad_getVmStatPages(vmStatText, 'Pages free');
inactivePages = matRad_getVmStatPages(vmStatText, 'Pages inactive');
speculativePages = matRad_getVmStatPages(vmStatText, 'Pages speculative');

if all(isnan([freePages inactivePages speculativePages]))
    return
end

availablePageValues = [freePages inactivePages speculativePages];
availablePages = sum(availablePageValues(~isnan(availablePageValues)));
availableBytes = availablePages * pageSize;

[status, totalText] = system('sysctl -n hw.memsize');
if status == 0
    totalBytes = str2double(strtrim(totalText));
end
end

function pages = matRad_getVmStatPages(vmStatText, label)
pages = NaN;
pattern = [regexptranslate('escape', label), ':\s+([0-9]+)\.'];
match = regexp(vmStatText, pattern, 'tokens', 'once');
if ~isempty(match)
    pages = str2double(match{1});
end
end

function [availableBytes, totalBytes, source] = matRad_getUnixSystemMemoryBytes()
[availableBytes, totalBytes, source] = matRad_getProcMemInfoMemoryBytes();
if ~isempty(availableBytes)
    return
end

[availableBytes, totalBytes, source] = matRad_getBsdVmStatMemoryBytes();
end

function [availableBytes, totalBytes, source] = matRad_getProcMemInfoMemoryBytes()
availableBytes = [];
totalBytes = [];
source = '';

try
    memInfoText = fileread('/proc/meminfo');
catch
    return
end

availableKb = matRad_getMemInfoKb(memInfoText, 'MemAvailable');
totalKb = matRad_getMemInfoKb(memInfoText, 'MemTotal');

if isnan(availableKb)
    freeKb = matRad_getMemInfoKb(memInfoText, 'MemFree');
    buffersKb = matRad_getMemInfoKb(memInfoText, 'Buffers');
    cachedKb = matRad_getMemInfoKb(memInfoText, 'Cached');
    reclaimableKb = matRad_getMemInfoKb(memInfoText, 'SReclaimable');
    shmemKb = matRad_getMemInfoKb(memInfoText, 'Shmem');
    availableComponents = [freeKb buffersKb cachedKb reclaimableKb];
    availableKb = sum(availableComponents(~isnan(availableComponents)));
    if ~isnan(shmemKb)
        availableKb = availableKb - shmemKb;
    end
end

if ~isnan(availableKb)
    availableBytes = max(0, availableKb) * 1024;
end
if ~isnan(totalKb)
    totalBytes = totalKb * 1024;
end
if ~isempty(availableBytes)
    source = '/proc/meminfo';
end
end

function value = matRad_getKeyValueNumber(text, label)
value = NaN;
pattern = [regexptranslate('escape', label), '=([0-9]+)'];
match = regexp(text, pattern, 'tokens', 'once');
if ~isempty(match)
    value = str2double(match{1});
end
end

function valueKb = matRad_getMemInfoKb(memInfoText, label)
valueKb = NaN;
pattern = [regexptranslate('escape', label), ':\s+([0-9]+)\s+kB'];
match = regexp(memInfoText, pattern, 'tokens', 'once');
if ~isempty(match)
    valueKb = str2double(match{1});
end
end

function [availableBytes, totalBytes, source] = matRad_getBsdVmStatMemoryBytes()
availableBytes = [];
totalBytes = [];
source = '';

pageSize = matRad_getSysctlValueBytes('hw.pagesize');
if isempty(pageSize)
    pageSize = 4096;
end

[status, vmStatText] = system('vmstat -s');
if status ~= 0
    return
end

freePages = matRad_getBsdVmStatPages(vmStatText, 'pages free');
inactivePages = matRad_getBsdVmStatPages(vmStatText, 'pages inactive');
cachedPages = matRad_getBsdVmStatPages(vmStatText, 'pages cached');

if all(isnan([freePages inactivePages cachedPages]))
    return
end

availablePageValues = [freePages inactivePages cachedPages];
availablePages = sum(availablePageValues(~isnan(availablePageValues)));
availableBytes = availablePages * pageSize;

totalBytes = matRad_getSysctlValueBytes('hw.memsize');
if isempty(totalBytes)
    totalBytes = matRad_getSysctlValueBytes('hw.physmem');
end
if isempty(totalBytes)
    totalBytes = matRad_getSysctlValueBytes('hw.physmem64');
end
if isempty(totalBytes)
    totalBytes = matRad_getSysctlValueBytes('hw.realmem');
end

source = 'vmstat/sysctl';
end

function pages = matRad_getBsdVmStatPages(vmStatText, label)
pages = NaN;
pattern = ['([0-9]+)\s+', regexptranslate('escape', label)];
match = regexp(vmStatText, pattern, 'tokens', 'once');
if ~isempty(match)
    pages = str2double(match{1});
end
end

function valueBytes = matRad_getSysctlValueBytes(key)
valueBytes = [];

[status, valueText] = system(['sysctl -n ', key]);
if status ~= 0
    return
end

valueBytes = str2double(strtrim(valueText));
if ~isfinite(valueBytes)
    valueBytes = [];
end
end
