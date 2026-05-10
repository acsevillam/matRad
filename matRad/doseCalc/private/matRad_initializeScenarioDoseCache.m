function cacheContext = matRad_initializeScenarioDoseCache(cfg,ctx,quantity,stf,matRad_cfg,calculationName,signatureExtras)
% matRad_initializeScenarioDoseCache creates a hashed streaming cache folder
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

if nargin < 6 || isempty(calculationName)
    calculationName = 'scenario dose';
end
if nargin < 7 || isempty(signatureExtras)
    signatureExtras = struct();
end

cacheRoot = cfg.CacheRoot;
if exist(cacheRoot,'file') == 2
    matRad_cfg.dispError('Could not create %s cache root "%s": a file already exists at that path.', ...
        calculationName,cacheRoot);
end
if exist(cacheRoot,'dir') ~= 7
    [success,message] = mkdir(cacheRoot);
    if ~success
        matRad_cfg.dispError('Could not create %s cache root "%s": %s', ...
            calculationName,cacheRoot,message);
    end
end

signature = buildCacheSignature(cfg,ctx,quantity,stf,calculationName, ...
    signatureExtras);
runHash = buildCacheRunHash(signature);
runDir = fullfile(cacheRoot,runHash);
counter = 0;
while exist(runDir,'dir') == 7
    counter = counter + 1;
    runDir = fullfile(cacheRoot,sprintf('%s_%03d',runHash,counter));
end

[success,message] = mkdir(runDir);
if ~success
    matRad_cfg.dispError('Could not create %s cache folder "%s": %s', ...
        calculationName,runDir,message);
end

cacheContext.root = cacheRoot;
cacheContext.runDir = runDir;
cacheContext.runHash = runHash;
cacheContext.created = true;
cacheContext.calculationName = calculationName;
cacheContext.signature = signature;
metadata = signature;
metadataFile = fullfile(runDir,'metadata.mat');
try
    save(metadataFile,'metadata','-v7');
catch ME
    cleanupCreatedRunDir(runDir);
    matRad_cfg.dispError('Could not write %s cache metadata "%s": %s', ...
        calculationName,metadataFile,ME.message);
end
end

function signature = buildCacheSignature(cfg,ctx,quantity,stf,calculationName,signatureExtras)
signature = struct();
signature.calculationName = calculationName;
if isfield(cfg,'IntervalMode')
    signature.intervalMode = cfg.IntervalMode;
end
if isfield(cfg,'RadiusMode')
    signature.radiusMode = cfg.RadiusMode;
end
signature.secondPassStrategy = cfg.SecondPassStrategy;
signature.quantity = quantity.name;
signature.quantityField = quantity.field;
signature.quantityScale = quantity.scale;
signature.refScen = cfg.refScen;
signature.scenarioDijIx = ctx.scenarioDijIx(:);
signature.scenarioCtScenIds = ctx.scenarioCtScenIds(:);
signature.scenarioWeights = ctx.scenarioWeights(:);
if isfield(ctx,'scenarioIds')
    signature.scenarioIds = ctx.scenarioIds(:);
end
signature.numVoxels = ctx.numVoxels;
signature.numBixels = ctx.numBixels;
signature.doseGrid = ctxDoseGridSignature(ctx);
signature.stf = stfSignature(stf);
signature.extras = signatureExtras;
signature.nonce = buildRunNonce();
end

function doseGrid = ctxDoseGridSignature(ctx)
doseGrid = struct();
doseGrid.numVoxels = ctx.numVoxels;
doseGrid.numBixels = ctx.numBixels;
if isfield(ctx,'doseGrid')
    doseGrid.doseGrid = ctx.doseGrid;
end
if isfield(ctx,'ctGrid')
    doseGrid.ctGrid = ctx.ctGrid;
end
end

function signature = stfSignature(stf)
signature = struct();
signature.numBeams = numel(stf);
if isempty(stf)
    signature.beams = [];
    return;
end

beams = struct('machine',{},'radiationMode',{},'numOfRays',{}, ...
    'totalNumOfBixels',{},'gantryAngle',{},'couchAngle',{},'isoCenter',{});
for i = 1:numel(stf)
    beams(i).machine = getStructFieldOrEmpty(stf(i),'machine');
    beams(i).radiationMode = getStructFieldOrEmpty(stf(i),'radiationMode');
    beams(i).numOfRays = getStructFieldOrEmpty(stf(i),'numOfRays');
    beams(i).totalNumOfBixels = getStructFieldOrEmpty(stf(i),'totalNumOfBixels');
    beams(i).gantryAngle = getStructFieldOrEmpty(stf(i),'gantryAngle');
    beams(i).couchAngle = getStructFieldOrEmpty(stf(i),'couchAngle');
    beams(i).isoCenter = getStructFieldOrEmpty(stf(i),'isoCenter');
end
signature.beams = beams;
end

function value = getStructFieldOrEmpty(s,fieldName)
if isfield(s,fieldName)
    value = s.(fieldName);
else
    value = [];
end
end

function runHash = buildCacheRunHash(signature)
text = serializeForHash(signature);
hashText = fnv1a32(text);
timestamp = datestr(now,'yyyymmddTHHMMSS');
runHash = lower(sprintf('%s_%s',timestamp,hashText));
end

function nonce = buildRunNonce()
clockVector = clock;
nonce = sprintf('%.17g_%04.0f%02.0f%02.0fT%02.0f%02.0f%012.6f_%.17g', ...
    now,clockVector(1),clockVector(2),clockVector(3),clockVector(4), ...
    clockVector(5),clockVector(6),cputime);
end

function cleanupCreatedRunDir(runDir)
if exist(runDir,'dir') == 7
    [~,~] = rmdir(runDir,'s');
end
end

function text = serializeForHash(value)
if isstruct(value)
    if numel(value) ~= 1
        parts = cell(numel(value),1);
        for i = 1:numel(value)
            parts{i} = serializeForHash(value(i));
        end
        text = ['[' strjoinCompat(parts,';') ']'];
    else
        fields = sort(fieldnames(value));
        parts = cell(numel(fields),1);
        for i = 1:numel(fields)
            parts{i} = [fields{i} '=' serializeForHash(value.(fields{i}))];
        end
        text = ['{' strjoinCompat(parts,';') '}'];
    end
elseif iscell(value)
    parts = cell(numel(value),1);
    for i = 1:numel(value)
        parts{i} = serializeForHash(value{i});
    end
    text = ['[' strjoinCompat(parts,';') ']'];
elseif isnumeric(value) || islogical(value)
    text = mat2str(value);
elseif ischar(value)
    text = value;
elseif isstring(value)
    text = char(value);
else
    text = class(value);
end
end

function hashText = fnv1a32(text)
bytes = uint8(text);
hashValue = uint32(2166136261);
prime = uint32(16777619);
for i = 1:numel(bytes)
    hashValue = bitxor(hashValue,uint32(bytes(i)));
    hashValue = uint32(mod(double(hashValue) * double(prime),2^32));
end
hashText = dec2hex(hashValue,8);
end

function out = strjoinCompat(parts,delimiter)
if isempty(parts)
    out = '';
    return;
end
out = parts{1};
for i = 2:numel(parts)
    out = [out delimiter parts{i}]; %#ok<AGROW>
end
end
