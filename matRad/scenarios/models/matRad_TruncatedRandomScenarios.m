classdef matRad_TruncatedRandomScenarios < matRad_RandomScenariosAbstract
    %  matRad_TruncatedRandomScenarios
    %  Implements randomly sampled scenarios truncated to the wcSigma
    %  normalized uncertainty radius equivalent to a diagonal Mahalanobis
    %  distance in the active uncertainty component space.
    %
    % constructor
    %   matRad_TruncatedRandomScenarios()
    %   matRad_TruncatedRandomScenarios(ct)
    %
    % input:
    %   ct:                 ct cube
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

    properties (SetAccess = protected)
        name = 'Truncated Randomly sampled Scenarios'
        shortName = 'truncatedRndScen'
    end

    methods

        function this = matRad_TruncatedRandomScenarios(ct)
            if nargin == 0
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end

            this@matRad_RandomScenariosAbstract(superclassArgs{:});
        end

    end

    methods (Access = protected)

        function scenarios = sampleScenarioValues(this, nScenarioSamples, components, activeIx)
            scenarios = zeros(nScenarioSamples, numel(components));
            if ~any(activeIx) || this.wcSigma <= 0
                return
            end

            numActiveComponents = sum(activeIx);
            randomValues = this.sampleTruncatedStandardNormalValues(nScenarioSamples, numActiveComponents);
            scales = [components(activeIx).scale];
            scenarios(:, activeIx) = bsxfun(@times, randomValues, scales);
        end

        function randomValues = sampleTruncatedStandardNormalValues(this, numRows, numColumns)
            if isempty(this.randomSeed)
                randomValues = this.drawTruncatedStandardNormalValues(numRows, numColumns);
            else
                rngState = rng;
                cleanupObj = onCleanup(@() rng(rngState));
                rng(this.randomSeed);
                randomValues = this.drawTruncatedStandardNormalValues(numRows, numColumns);
                clear cleanupObj;
            end
        end

        function randomValues = drawTruncatedStandardNormalValues(this, numRows, numColumns)
            randomValues = zeros(numRows, numColumns);
            numAccepted = 0;
            numBatches = 0;
            maxNumBatches = 1000;
            batchSize = max(100, numRows);
            tolerance = 100 * eps(max(1, this.wcSigma));

            while numAccepted < numRows && numBatches < maxNumBatches
                candidates = randn(batchSize, numColumns);
                normalizedRadius = sqrt(sum(candidates.^2, 2));
                accepted = candidates(normalizedRadius <= this.wcSigma + tolerance, :);
                numNew = min(size(accepted, 1), numRows - numAccepted);

                if numNew > 0
                    targetIx = (numAccepted + 1):(numAccepted + numNew);
                    randomValues(targetIx, :) = accepted(1:numNew, :);
                    numAccepted = numAccepted + numNew;
                end

                numBatches = numBatches + 1;
            end

            if numAccepted < numRows
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError(['Could not sample %d truncated random scenarios within wcSigma %.3g. ' ...
                                      'Increase wcSigma or reduce the number of active uncertainty components.'], ...
                                     numRows, this.wcSigma);
            end
        end

    end

end
