classdef matRad_RegressionProbeEngine < DoseEngines.matRad_DoseEngineBase
% Minimal concrete engine exposing base-engine initialization diagnostics.

    properties (Constant)
        shortName = 'regressionProbe';
        name = 'Regression probe';
        possibleRadiationModes = {'photons'};
    end

    methods
        function this = matRad_RegressionProbeEngine(pln)
            this@DoseEngines.matRad_DoseEngineBase(pln);
        end

        function applicators = supportedScenarioApplicators(~)
            applicators = {'ct','setup','range'};
        end
    end

    methods (Access = protected)
        function dij = calcDose(this,ct,cst,stf)
            dij = this.initDoseCalc(ct,cst,stf);
            dij.regressionProbe.numOfColumnsDij = this.numOfColumnsDij;
            dij.regressionProbe.VctGrid = this.VctGrid;
            dij.regressionProbe.VdoseGrid = this.VdoseGrid;
            dij.regressionProbe.VctGridScenIx = this.VctGridScenIx;
            dij.regressionProbe.VdoseGridScenIx = this.VdoseGridScenIx;
            dij.regressionProbe.VctGridMask = this.VctGridMask;
            dij.regressionProbe.VdoseGridMask = this.VdoseGridMask;
            dij.regressionProbe.robustVoxelsOnGrid = this.robustVoxelsOnGrid;
            dij = this.finalizeDose(dij);
        end
    end

    methods (Static)
        function [available,msg] = isAvailable(~,~)
            available = true;
            msg = '';
        end
    end
end
