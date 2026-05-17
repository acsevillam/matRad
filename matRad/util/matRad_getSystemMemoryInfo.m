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
p.parse(varargin{:});

memoryInfo = struct('availableBytes', [], 'totalBytes', [], ...
                    'reserveBytes', [], 'usableBytes', [], 'source', '');

[availableBytes, totalBytes, source] = matRad_getAvailableSystemMemoryBytes();

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
