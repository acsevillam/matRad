function optimizer = matRad_createOptimizer(optimizerName,varargin)
% matRad_createOptimizer creates a configured optimizer by name.
%
% call
%   optimizer = matRad_createOptimizer(optimizerName)
%   optimizer = matRad_createOptimizer(optimizerName,'fallbackOptimizer','IPOPT')
%
% input
%   optimizerName:      Name of the optimizer implementation.
%   varargin:           optional Name-Value pair arguments
%
% input (optional Name-Value pairs)
%   fallbackOptimizer:  Optimizer used when optimizerName is unknown. Empty
%                       keeps strict error behavior.
%
% output
%   optimizer:          Optimizer instance.
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

matRad_cfg = MatRad_Config.instance();
options = parseOptions(varargin{:});

if nargin < 1 || isempty(optimizerName)
    optimizerName = 'IPOPT';
end
optimizerName = normalizeOptimizerName(optimizerName,'optimizerName',matRad_cfg);

[optimizerConstructor,optimizerAvailabilityCheck,optimizerName] = ...
    resolveOptimizer(optimizerName);
if isempty(optimizerConstructor)
    if isempty(options.fallbackOptimizer)
        matRad_cfg.dispError('Optimizer ''%s'' not known.',optimizerName);
    end
    unknownOptimizerName = optimizerName;
    [optimizerConstructor,optimizerAvailabilityCheck,optimizerName] = ...
        resolveOptimizer(options.fallbackOptimizer);
    if isempty(optimizerConstructor)
        matRad_cfg.dispError('Fallback optimizer ''%s'' not known.', ...
            options.fallbackOptimizer);
    end
    matRad_cfg.dispWarning('Optimizer ''%s'' not known! Fallback to %s!', ...
        unknownOptimizerName,optimizerName);
end

if ~optimizerAvailabilityCheck()
    matRad_cfg.dispError('Optimizer ''%s'' not available.',optimizerName);
end

optimizer = optimizerConstructor();
end

function options = parseOptions(varargin)
matRad_cfg = MatRad_Config.instance();
options = struct('fallbackOptimizer','');

if mod(numel(varargin),2) ~= 0
    matRad_cfg.dispError('Optimizer options must be name-value pairs.');
end

for i = 1:2:numel(varargin)
    name = normalizeOptionName(varargin{i},matRad_cfg);
    value = varargin{i+1};
    switch lower(name)
        case 'fallbackoptimizer'
            if isempty(value)
                options.fallbackOptimizer = '';
            else
                options.fallbackOptimizer = normalizeOptimizerName( ...
                    value,'fallbackOptimizer',matRad_cfg);
            end
        otherwise
            matRad_cfg.dispError('Unknown optimizer factory option ''%s''.',name);
    end
end
end

function name = normalizeOptionName(name,matRad_cfg)
if isstring(name) && isscalar(name)
    name = char(name);
end
if ~ischar(name)
    matRad_cfg.dispError('Optimizer option names must be character vectors.');
end
end

function optimizerName = normalizeOptimizerName(optimizerName,valueName,matRad_cfg)
if isstring(optimizerName) && isscalar(optimizerName)
    optimizerName = char(optimizerName);
end
if ~ischar(optimizerName)
    matRad_cfg.dispError('%s must be a character vector or string scalar.', ...
        valueName);
end
end

function [optimizerConstructor,optimizerAvailabilityCheck,optimizerName] = resolveOptimizer(optimizerName)
knownOptimizerNames = {'IPOPT','fmincon','simulannealbnd'};
knownOptimizerConstructors = {@matRad_OptimizerIPOPT, ...
    @matRad_OptimizerFmincon,@matRad_OptimizerSimulannealbnd};
knownOptimizerAvailabilityChecks = {@() matRad_OptimizerIPOPT.IsAvailable(), ...
    @() matRad_OptimizerFmincon.IsAvailable(), ...
    @() matRad_OptimizerSimulannealbnd.IsAvailable()};

optimizerConstructor = [];
optimizerAvailabilityCheck = [];

matchIx = find(strcmp(optimizerName,knownOptimizerNames),1);
if ~isempty(matchIx)
    optimizerName = knownOptimizerNames{matchIx};
    optimizerConstructor = knownOptimizerConstructors{matchIx};
    optimizerAvailabilityCheck = knownOptimizerAvailabilityChecks{matchIx};
end
end
