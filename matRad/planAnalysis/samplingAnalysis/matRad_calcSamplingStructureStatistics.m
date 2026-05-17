function cstStat = matRad_calcSamplingStructureStatistics(cst, caSampRes, scenWeights, percentiles)
% matRad_calcSamplingStructureStatistics calculates structure sampling stats
%
% call
%   cstStat = matRad_calcSamplingStructureStatistics(cst,caSampRes,scenWeights,percentiles)
%
% input
%   cst:         matRad cst cell array represented on the analysis CT grid
%   caSampRes:   cell array of sampling result structs
%   scenWeights: scenario weights for weighted statistics
%   percentiles: percentiles evaluated for DVH and QI statistics
%
% output
%   cstStat:     structure-wise DVH and QI sampling statistics
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

percentileNames = cell(numel(percentiles), 1);
for i = 1:numel(percentiles)
    percentileNames{i} = ['P', num2str(percentiles(i) * 100)];
end
metric = vertcat({'mean'; 'min'; 'max'; 'std'}, percentileNames{:});

cstStat = struct();
for i = 1:size(cst, 1)
    cstStat(i).name = caSampRes(1).qi(i).name;
    for l = 1:numel(caSampRes)
        if any(~strcmp(cstStat(i).name, {caSampRes(l).dvh(i).name, caSampRes(l).qi(i).name}))
            matRad_cfg.dispError('matRad: Error, wrong structure.');
        end
        cstStat(i).dvh(l).doseGrid = caSampRes(l).dvh(i).doseGrid;
        cstStat(i).dvh(l).volumePoints = caSampRes(l).dvh(i).volumePoints;
        cstStat(i).qi(l) = caSampRes(l).qi(i);
        cstStat(i).w(l) = scenWeights(l)';
    end

    cstStat(i).percentiles = percentiles;
    cstStat(i).metric = metric;
    cstStat(i).dvhStat = matRad_calcDvhStat(cstStat(i).dvh, cstStat(i).percentiles, cstStat(i).w);
    cstStat(i).qiStat = matRad_calcQiStat(cstStat(i).qi, cstStat(i).percentiles, cstStat(i).w, metric);
end

end

function dvhStat = matRad_calcDvhStat(dvh, percentiles, weights)
doseGrid = dvh(1).doseGrid;
dvhMat = NaN * ones(numel(dvh), numel(dvh(1).volumePoints));
for j = 1:numel(dvh)
    dvhMat(j, :) = dvh(j).volumePoints;
end
dvhMat(isnan(dvhMat)) = 0;

dvhStat.mean.doseGrid = doseGrid;
dvhStat.mean.volumePoints = matRad_weightedMean(dvhMat, weights);

dvhStat.min.doseGrid = doseGrid;
dvhStat.min.volumePoints = min(dvhMat);

dvhStat.max.doseGrid = doseGrid;
dvhStat.max.volumePoints = max(dvhMat);

dvhStat.std.doseGrid = doseGrid;
dvhStat.std.volumePoints = matRad_weightedStd(dvhMat, weights);

dvhStat.percDVH = NaN * ones(numel(percentiles), numel(doseGrid));
for j = 1:size(dvhMat, 2)
    dvhStat.percDVH(:, j) = matRad_weightedQuantile(dvhMat(:, j), percentiles, weights', false);
end
end

function qiStat = matRad_calcQiStat(qi, percentiles, weights, metric)
fields = fieldnames(qi);
if any(strcmp('VOIname', fields))
    qi = rmfield(qi, 'VOIname');
end
fields = fieldnames(qi);
qiStruct = qi;

qiStatH = struct();
for j = 1:numel(fields)
    if numel([qiStruct(:).(fields{j})]) == numel(weights)
        qiStatH(1).(fields{j}) = matRad_weightedMean([qiStruct(:).(fields{j})], weights);
        qiStatH(2).(fields{j}) = min([qiStruct(:).(fields{j})]);
        qiStatH(3).(fields{j}) = max([qiStruct(:).(fields{j})]);
        qiStatH(4).(fields{j}) = matRad_weightedStd([qiStruct(:).(fields{j})], weights);
        wQ = matRad_weightedQuantile([qiStruct(:).(fields{j})], percentiles, weights', false);
        for k = 1:numel(wQ)
            qiStatH(k + 4).(fields{j}) = wQ(k);
        end
    else
        for k = 1:(4 + numel(percentiles))
            qiStatH(k).(fields{j}) = [];
        end
    end
end

env = matRad_getEnvironment();
if strcmp(env, 'MATLAB')
    qiStat = struct2table(qiStatH);
    qiStat.Properties.RowNames = metric;
else
    qiStat = qiStatH;
end
end

function meanValue = matRad_weightedMean(values, weights)
if exist('weights', 'var') && ~isempty(weights)
    if isvector(values) && isvector(weights)
        meanValue = reshape(weights, 1, []) * reshape(values, [], 1) / sum(weights);
    else
        meanValue = reshape(weights, 1, []) * values ./ sum(weights);
    end
else
    meanValue = mean(values);
end
end

function stdValue = matRad_weightedStd(values, weights)
if exist('weights', 'var') && ~isempty(weights)
    if isvector(values) && isvector(weights)
        values = reshape(values, [], 1);
        weights = reshape(weights, [], 1);
        mu = sum(weights .* values) ./ sum(weights);
        stdValue = sqrt(sum(weights .* (values - mu).^2) ./ sum(weights));
    else
        mu = matRad_weightedMean(values, weights);
        stdValue = sqrt(reshape(weights, 1, []) * ...
                        (bsxfun(@minus, values, mu).^2) ./ sum(weights));
    end
else
    stdValue = std(values, 1);
end
end
