function matRad_showDVHFromSampling(caSamp,scale,cst,pln,scenarios,doseWindow,dvhType,refScen,lineStyleIndicator,varargin)
% matRad_showDVHFromSampling displays sampled dose-volume histograms
%
% call
%   matRad_showDVHFromSampling(caSamp,scale,cst,pln,scenarios)
%   matRad_showDVHFromSampling(caSamp,scale,cst,pln,scenarios,doseWindow,dvhType,refScen,lineStyleIndicator)
%   matRad_showDVHFromSampling(...,'ctScenProb',ctScenProb)
%   matRad_showDVHFromSampling(...,'scenWeights',scenWeights)
%
% input
%   caSamp:             sampled scenario result struct with per-scenario DVHs
%   scale:              dose scaling factor applied to the DVH dose grid
%   cst:                matRad cst struct
%   pln:                matRad pln struct with multScen scenario model
%   scenarios:          sampled scenario indices to plot
%   doseWindow:         (optional) 1x2 dose axis limits
%   dvhType:            (optional) 'multiscenario', 'minmax', or 'trustband'
%   refScen:            (optional) scenario index highlighted as reference
%   lineStyleIndicator: (optional) integer selecting the line style
%   ctScenProb:         (optional) CT scenario probability override as
%                       Name-Value pair [ctScenIndex probability]
%   scenWeights:        (optional) resolved scenario weights as Name-Value
%                       pair. If provided, these weights are used directly.
%
% output
%   graphical display of sampled DVH curves or bands
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

if ~exist('dvhType','var') || isempty(dvhType)
    dvhType = 'multiscenario';
end

if isstring(dvhType) && isscalar(dvhType)
    dvhType = char(dvhType);
end

validDvhTypes = {'multiscenario','minmax','trustband'};
if ~ischar(dvhType) || ~any(strcmp(dvhType,validDvhTypes))
    matRad_cfg.dispError('Unsupported dvhType. Use ''multiscenario'', ''minmax'', or ''trustband''.');
end

if ~exist('scale','var') || isempty(scale)
    scale = 1.0;
end

if ~exist('refScen','var') || isempty(refScen)
    refScen = 1;
end

if ~exist('lineStyleIndicator','var') || isempty(lineStyleIndicator)
    lineStyleIndicator = 1;
end

p = inputParser;
p.CaseSensitive = false;
p.addParameter('ctScenProb',[],@(x) isempty(x) || (isnumeric(x) && ismatrix(x) && size(x,2) == 2 && all(isfinite(x(:))) && all(x(:,2) >= 0)));
p.addParameter('scenWeights',[],@(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.parse(varargin{:});
ctScenProb = p.Results.ctScenProb;
scenWeights = p.Results.scenWeights;

if ~isempty(ctScenProb) && ~isempty(scenWeights)
    matRad_cfg.dispError('Specify either ctScenProb or scenWeights, not both.');
end

% create new figure and set default line style indicator if not explictly
% specified
hold on;

cstNames = cst(:,2);
cstInfo = cst(:,5);

numOfVois = numel(cstNames);

%% print the dvh
try
    colorMx = cellfun(@(c) c.visibleColor,cstInfo,'UniformOutput',false);
    colorMx = cell2mat(colorMx);
catch
    colorMx    = colorcube;
    colorMx    = colorMx(1:floor(64/numOfVois):64,:);
end

lineStyles = {'-',':','--','-.'};

maxDVHvol  = 0;
maxDVHdose = 0;

switch dvhType

    case 'multiscenario'

        for i = 1:numOfVois
            for k = scenarios
                if cst{i,5}.Visible
                    % cut off at the first zero value where there is no more signal
                    % behind
                    ix      = lastDvhPlotIndex(caSamp(k).dvh(i).volumePoints);
                    currDvh = [caSamp(k).dvh(i).doseGrid(1:ix)*scale;caSamp(k).dvh(i).volumePoints(1:ix)];

                    p=plot(currDvh(1,:),currDvh(2,:),'LineWidth',0.5,'Color',[colorMx(i,:) 0.4], ...
                        'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cst{i,2});

                    if(k==refScen)
                        p.LineWidth = 2;
                        p.Color = [colorMx(i,:) 1];
                        p.Annotation.LegendInformation.IconDisplayStyle = 'on';
                    else
                        p.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    end

                    maxDVHvol  = max(maxDVHvol,max(currDvh(2,:)));
                    maxDVHdose = max(maxDVHdose,max(currDvh(1,:)));

                end
            end
        end

    case 'minmax'
        for i = 1:numOfVois
            if cst{i,5}.Visible
                doseGrid = caSamp(1).dvh(i).doseGrid*scale;
                allDvh = zeros(numel(scenarios) + 1,numel(doseGrid));
                allDvh(1,:) = doseGrid;

                for s = 1:numel(scenarios)
                    k = scenarios(s);

                    if(k==refScen)
                        % cut off at the first zero value where there is no more signal
                        % behind
                        ix      = lastDvhPlotIndex(caSamp(k).dvh(i).volumePoints);
                        refDvh = [caSamp(k).dvh(i).doseGrid(1:ix)*scale;caSamp(k).dvh(i).volumePoints(1:ix)];

                        p=plot(refDvh(1,:),refDvh(2,:),'LineWidth',0.5,'Color',[colorMx(i,:) 0.4], ...
                            'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cst{i,2});

                        p.LineWidth = 2;
                        p.Color = [colorMx(i,:) 1];
                        p.Annotation.LegendInformation.IconDisplayStyle = 'on';
                    end

                    allDvh(s + 1,:) = caSamp(k).dvh(i).volumePoints;

                end

                bandDvh = [allDvh(1,:);min(allDvh(2:end,:));max(allDvh(2:end,:))];

                f = fill([bandDvh(1,:) flip(bandDvh(1,:))],[bandDvh(2,:) flip(bandDvh(3,:))],colorMx(i,:));
                f.FaceAlpha = 0.2;
                f.FaceColor = colorMx(i,:);
                f.EdgeColor = colorMx(i,:);
                f.LineWidth = 0.05;
                f.Annotation.LegendInformation.IconDisplayStyle = 'off';

                maxDVHvol  = max(maxDVHvol,max(bandDvh(3,:)));
                maxDVHdose = max(maxDVHdose,max(bandDvh(1,:)));

            end
        end

    case 'trustband'

        scenWeights = matRad_getSamplingScenarioWeights(pln,numel(caSamp),ctScenProb,scenWeights);

        for i = 1:numOfVois

            if cst{i,5}.Visible

                doseGrid = caSamp(1).dvh(i).doseGrid*scale;
                allDvh = zeros(numel(scenarios) + 1,numel(doseGrid));
                allDvh(1,:) = doseGrid;

                for s = 1:numel(scenarios)
                    k = scenarios(s);
                    if(k==refScen)
                        % cut off at the first zero value where there is no more signal
                        % behind
                        ix      = lastDvhPlotIndex(caSamp(k).dvh(i).volumePoints);
                        refDvh = [caSamp(k).dvh(i).doseGrid(1:ix)*scale;caSamp(k).dvh(i).volumePoints(1:ix)];

                        p=plot(refDvh(1,:),refDvh(2,:),'LineWidth',0.5,'Color',[colorMx(i,:) 0.4], ...
                            'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cst{i,2});

                        p.LineWidth = 2;
                        p.Color = [colorMx(i,:) 1];
                        p.Annotation.LegendInformation.IconDisplayStyle = 'on';
                    end

                    allDvh(s + 1,:) = caSamp(k).dvh(i).volumePoints;
                end

                currScenWeights = getCurrentScenarioWeights(scenWeights,scenarios,numel(caSamp));
                dvhValues = allDvh(2:end,:);
                meanDVHVolW = wMean(dvhValues,currScenWeights);
                stdDVHVolW  = wStd(dvhValues,currScenWeights);

                minDvh=max(meanDVHVolW-stdDVHVolW,0);
                maxDvh=min(meanDVHVolW+stdDVHVolW,100);

                bandDvh = [allDvh(1,:);minDvh;maxDvh];

                p=plot(bandDvh(1,:),meanDVHVolW,'LineWidth',0.5,'Color',[colorMx(i,:) 0.4], ...
                    'LineStyle',lineStyles{2},'DisplayName',cst{i,2});

                p.LineWidth = 2;
                p.Color = [colorMx(i,:) 1];
                p.Annotation.LegendInformation.IconDisplayStyle = 'off';

                f = fill([bandDvh(1,:) flip(bandDvh(1,:))],[bandDvh(2,:) flip(bandDvh(3,:))],colorMx(i,:));
                f.FaceAlpha = 0.2;
                f.FaceColor = colorMx(i,:);
                f.EdgeColor = colorMx(i,:);
                f.LineWidth = 0.05;
                f.Annotation.LegendInformation.IconDisplayStyle = 'off';

                maxDVHvol  = max(maxDVHvol,max([max(bandDvh(2,:)) max(bandDvh(3,:))]));
                maxDVHdose = max(maxDVHdose,max(bandDvh(1,:)));

            end
        end

end

fontSizeValue = 14;
myLegend = legend('show','location','NorthEast');
set(myLegend,'FontSize',10,'Interpreter','none');
legend boxoff

if ~exist('doseWindow', 'var') || isempty(doseWindow)
    doseWindow = [0 1.4*maxDVHdose];
else
    doseWindow = doseWindow(:)';
end

if numel(doseWindow) < 2 || ~all(isfinite(doseWindow(1:2))) || doseWindow(2) <= doseWindow(1)
    fallbackMaxDose = 1.4*maxDVHdose;
    if ~isfinite(fallbackMaxDose) || fallbackMaxDose <= 0
        fallbackMaxDose = 1;
    end
    doseWindow = [0 fallbackMaxDose];
end

if maxDVHvol <= 0
    maxDVHvol = 1;
end

ylim([0 1.05*maxDVHvol]);
xlim(doseWindow(1:2));

grid on,grid minor
box(gca,'on');
set(gca,'LineWidth',.5,'FontSize',fontSizeValue);
ylabel('Volume [%]','FontSize',fontSizeValue)

if strcmp(pln.bioParam.model,'none')
    xlabel('Dose [Gy]','FontSize',fontSizeValue);
else
    xlabel('RBE x Dose [Gy(RBE)]','FontSize',fontSizeValue);
end

end

function currScenWeights = getCurrentScenarioWeights(scenWeights,scenarios,numScenarios)
    matRad_cfg = MatRad_Config.instance();

    if numel(scenWeights) == numScenarios
        currScenWeights = scenWeights(scenarios(:));
    elseif numel(scenWeights) == numel(scenarios)
        currScenWeights = scenWeights(:);
    else
        matRad_cfg.dispError('Number of scenario weights does not match selected scenarios.');
    end
end

function ix = lastDvhPlotIndex(volumePoints)
    lastPositiveIx = find(volumePoints > 0,1,'last');
    if isempty(lastPositiveIx)
        ix = 1;
    else
        ix = min(lastPositiveIx + 1,numel(volumePoints));
    end
end

function meanValue = wMean(values,weights)
    weights = weights(:);
    meanValue = weights' * values ./ sum(weights);
end

function stdValue = wStd(values,weights)
    meanValue = wMean(values,weights);
    stdValue = sqrt(wMean(bsxfun(@minus,values,meanValue).^2,weights));
end
