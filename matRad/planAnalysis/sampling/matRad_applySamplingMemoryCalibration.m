function [samplingMemoryEstimate, calibration] = matRad_applySamplingMemoryCalibration(samplingMemoryEstimate, calibration, samplingResourceConfig)
% matRad_applySamplingMemoryCalibration applies safe sampling memory calibration
%
% call:
%   [samplingMemoryEstimate,calibration] =
%       matRad_applySamplingMemoryCalibration(samplingMemoryEstimate,calibration,samplingResourceConfig)
%
% input:
%   samplingMemoryEstimate: sampling memory estimate struct
%   calibration:            sampling calibration measurement struct
%   samplingResourceConfig: sampling resource configuration struct
%
% output:
%   samplingMemoryEstimate: sampling memory estimate after calibration policy
%   calibration:            calibration metadata with selected action
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calibration.allowReduction = ...
    samplingResourceConfig.allowCalibrationToReduceWorkerMemory;
calibration.calibratedMinWorkerBytes = ...
    samplingResourceConfig.calibratedMinForwardDoseWorkerMemoryBytes;
calibration.minReliableMeasuredBytes = ...
    samplingResourceConfig.calibrationMinReliableMeasuredBytes;
calibration.minReliableReductionRatio = ...
    samplingResourceConfig.calibrationMinReliableReductionRatio;

if ~calibration.measurementReliable
    calibration.action = 'unreliable';
    return
end

staticWorkerBytes = samplingMemoryEstimate.rawWorkerBytes;
calibratedWorkerBytes = max(calibration.measuredWorkerBytes, ...
    samplingResourceConfig.calibratedMinForwardDoseWorkerMemoryBytes);
calibration.calibratedWorkerBytes = calibratedWorkerBytes;
calibration.action = 'kept';

if calibratedWorkerBytes > staticWorkerBytes
    samplingMemoryEstimate.rawWorkerBytes = calibratedWorkerBytes;
    samplingMemoryEstimate.estimateBasis = 'samplingForwardDoseCalibrated';
    calibration.usedForPlanning = true;
    calibration.action = 'raised';
elseif samplingResourceConfig.allowCalibrationToReduceWorkerMemory && ...
        calibratedWorkerBytes < staticWorkerBytes
    if matRad_samplingCalibrationReductionIsReliable( ...
            calibration.measuredWorkerBytes, calibratedWorkerBytes, ...
            staticWorkerBytes, samplingResourceConfig)
        samplingMemoryEstimate.rawWorkerBytes = calibratedWorkerBytes;
        samplingMemoryEstimate.estimateBasis = 'samplingForwardDoseCalibrated';
        calibration.usedForPlanning = true;
        calibration.action = 'lowered';
    else
        calibration.action = 'undermeasured';
        calibration.undermeasurementReason = 'calibration_delta_too_small';
    end
end

end

function isReliable = matRad_samplingCalibrationReductionIsReliable(measuredWorkerBytes, ...
                                                                    calibratedWorkerBytes, ...
                                                                    staticWorkerBytes, ...
                                                                    samplingResourceConfig)
isReliable = measuredWorkerBytes >= ...
    samplingResourceConfig.calibrationMinReliableMeasuredBytes && ...
    calibratedWorkerBytes >= ...
    samplingResourceConfig.calibrationMinReliableReductionRatio * staticWorkerBytes;
end
