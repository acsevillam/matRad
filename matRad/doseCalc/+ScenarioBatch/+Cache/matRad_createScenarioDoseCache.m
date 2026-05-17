function cacheContext = matRad_createScenarioDoseCache(cfg, ctx, quantity, stf, matRadCfg, calculationName, signatureExtras)
% matRad_createScenarioDoseCache creates a hashed scenario-batch cache folder
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
if exist(cacheRoot, 'file') == 2
    matRadCfg.dispError('Could not create %s cache root "%s": a file already exists at that path.', ...
                        calculationName, cacheRoot);
end
if exist(cacheRoot, 'dir') ~= 7
    [success, message] = mkdir(cacheRoot);
    if ~success
        matRadCfg.dispError('Could not create %s cache root "%s": %s', ...
                            calculationName, cacheRoot, message);
    end
end

signature = matRad_buildCacheSignature(cfg, ctx, quantity, stf, calculationName, ...
                                       signatureExtras);
runHash = matRad_buildCacheRunHash(signature);
runDir = fullfile(cacheRoot, runHash);
counter = 0;
while exist(runDir, 'dir') == 7
    counter = counter + 1;
    runDir = fullfile(cacheRoot, sprintf('%s_%03d', runHash, counter));
end

[success, message] = mkdir(runDir);
if ~success
    matRadCfg.dispError('Could not create %s cache folder "%s": %s', ...
                        calculationName, runDir, message);
end

cacheContext.root = cacheRoot;
cacheContext.runDir = runDir;
cacheContext.runHash = runHash;
cacheContext.created = true;
cacheContext.calculationName = calculationName;
cacheContext.signature = signature;
metadata = signature;
metadataFile = fullfile(runDir, 'metadata.mat');
try
    save(metadataFile, 'metadata', '-v7');
catch ME
    matRad_cleanupCreatedRunDir(runDir);
    matRadCfg.dispError('Could not write %s cache metadata "%s": %s', ...
                        calculationName, metadataFile, ME.message);
end
end

function signature = matRad_buildCacheSignature(cfg, ctx, quantity, stf, calculationName, signatureExtras)
signature = struct();
signature.calculationName = calculationName;
signature.secondPassStrategy = cfg.SecondPassStrategy;
signature.quantity = quantity.name;
signature.quantityField = quantity.field;
signature.quantityScale = quantity.scale;
signature.refScen = cfg.refScen;
signature.scenarioDijIx = ctx.scenarioDijIx(:);
signature.scenarioCtScenIds = ctx.scenarioCtScenIds(:);
signature.scenarioWeights = ctx.scenarioWeights(:);
if isfield(ctx, 'scenarioIds')
    signature.scenarioIds = ctx.scenarioIds(:);
end
signature.numVoxels = ctx.numVoxels;
signature.numBixels = ctx.numBixels;
signature.doseGrid = matRad_ctxDoseGridSignature(ctx);
signature.stf = matRad_stfSignature(stf);
signature = matRad_promoteSignatureExtras(signature, signatureExtras);
signature.extras = signatureExtras;
signature.nonce = matRad_buildRunNonce();
end

function signature = matRad_promoteSignatureExtras(signature, signatureExtras)
extraNames = fieldnames(signatureExtras);
for i = 1:numel(extraNames)
    fieldName = extraNames{i};
    if ~isfield(signature, fieldName)
        signature.(fieldName) = signatureExtras.(fieldName);
    end
end
end

function doseGrid = matRad_ctxDoseGridSignature(ctx)
doseGrid = struct();
doseGrid.numVoxels = ctx.numVoxels;
doseGrid.numBixels = ctx.numBixels;
if isfield(ctx, 'doseGrid')
    doseGrid.doseGrid = ctx.doseGrid;
end
if isfield(ctx, 'ctGrid')
    doseGrid.ctGrid = ctx.ctGrid;
end
end

function signature = matRad_stfSignature(stf)
signature = struct();
signature.numBeams = numel(stf);
if isempty(stf)
    signature.beams = [];
    return
end

beams = struct('machine', {}, 'radiationMode', {}, 'numOfRays', {}, ...
               'totalNumOfBixels', {}, 'gantryAngle', {}, 'couchAngle', {}, 'isoCenter', {});
for i = 1:numel(stf)
    beams(i).machine = matRad_getStructFieldOrEmpty(stf(i), 'machine');
    beams(i).radiationMode = matRad_getStructFieldOrEmpty(stf(i), 'radiationMode');
    beams(i).numOfRays = matRad_getStructFieldOrEmpty(stf(i), 'numOfRays');
    beams(i).totalNumOfBixels = matRad_getStructFieldOrEmpty(stf(i), 'totalNumOfBixels');
    beams(i).gantryAngle = matRad_getStructFieldOrEmpty(stf(i), 'gantryAngle');
    beams(i).couchAngle = matRad_getStructFieldOrEmpty(stf(i), 'couchAngle');
    beams(i).isoCenter = matRad_getStructFieldOrEmpty(stf(i), 'isoCenter');
end
signature.beams = beams;
end

function value = matRad_getStructFieldOrEmpty(s, fieldName)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = [];
end
end

function runHash = matRad_buildCacheRunHash(signature)
text = matRad_serializeForHash(signature);
hashText = matRad_fnv1a32(text);
timestamp = datestr(now, 'yyyymmddTHHMMSS');
runHash = lower(sprintf('%s_%s', timestamp, hashText));
end

function nonce = matRad_buildRunNonce()
clockVector = clock;
nonce = sprintf('%.17g_%04.0f%02.0f%02.0fT%02.0f%02.0f%012.6f_%.17g', ...
                now, clockVector(1), clockVector(2), clockVector(3), clockVector(4), ...
                clockVector(5), clockVector(6), cputime);
end

function matRad_cleanupCreatedRunDir(runDir)
if exist(runDir, 'dir') == 7
    [~, ~] = rmdir(runDir, 's');
end
end

function text = matRad_serializeForHash(value)
if isstruct(value)
    if numel(value) ~= 1
        parts = cell(numel(value), 1);
        for i = 1:numel(value)
            parts{i} = matRad_serializeForHash(value(i));
        end
        text = ['[' matRad_strjoinCompat(parts, ';') ']'];
    else
        fields = sort(fieldnames(value));
        parts = cell(numel(fields), 1);
        for i = 1:numel(fields)
            parts{i} = [fields{i} '=' matRad_serializeForHash(value.(fields{i}))];
        end
        text = ['{' matRad_strjoinCompat(parts, ';') '}'];
    end
elseif iscell(value)
    parts = cell(numel(value), 1);
    for i = 1:numel(value)
        parts{i} = matRad_serializeForHash(value{i});
    end
    text = ['[' matRad_strjoinCompat(parts, ';') ']'];
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

function hashText = matRad_fnv1a32(text)
bytes = uint8(text);
hashValue = uint32(2166136261);
prime = uint32(16777619);
for i = 1:numel(bytes)
    hashValue = bitxor(hashValue, uint32(bytes(i)));
    hashValue = uint32(mod(double(hashValue) * double(prime), 2^32));
end
hashText = dec2hex(hashValue, 8);
end

function out = matRad_strjoinCompat(parts, delimiter)
if isempty(parts)
    out = '';
    return
end
out = parts{1};
for i = 2:numel(parts)
    out = [out delimiter parts{i}]; %#ok<AGROW>
end
end
