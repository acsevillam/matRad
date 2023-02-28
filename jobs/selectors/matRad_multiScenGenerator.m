function [multScen] = matRad_multiScenGenerator(scen_mode,run_config,mode,ct)

if(mode=="optimization")
    wcFactor = run_config.wcFactor;
end

if(mode=="sampling")
    wcFactor = run_config.sampling_wcFactor;
end

shiftSD=run_config.shiftSD;

switch scen_mode
    case "nomScen"
        multScen = matRad_multScen(ct,'nomScen');
    case "wcScen"
        multScen = matRad_multScen(ct,'wcScen');
        multScen.wcFactor=wcFactor;
        multScen.numOfShiftScen = [3 3 3];
        multScen.shiftSD = shiftSD;
        multScen.includeNomScen=true;
        multScen.numOfRangeShiftScen=6;
        multScen.rangeRelSD=0;
        multScen.rangeAbsSD=0;
        multScen.scenCombType = 'combined';
    case "impScen5"
        multScen = matRad_multScen(ct,'impScen');
        multScen.wcFactor=wcFactor;
        multScen.numOfShiftScen = [5 5 5];
        multScen.shiftSD = shiftSD;
        multScen.shiftGenType = 'equidistant';
        multScen.shiftCombType='individual';
        multScen.numOfRangeShiftScen=12;
        multScen.rangeRelSD=0;
        multScen.rangeAbsSD=0;
        multScen.scenCombType = 'combined';
        multScen.includeNomScen=true;
    case "impScen7"
        multScen = matRad_multScen(ct,'impScen');
        multScen.wcFactor=wcFactor;
        multScen.numOfShiftScen = [7 7 7];
        multScen.shiftSD = shiftSD;
        multScen.shiftGenType = 'equidistant';
        multScen.shiftCombType='individual';
        multScen.numOfRangeShiftScen=12;
        multScen.rangeRelSD=0;
        multScen.rangeAbsSD=0;
        multScen.scenCombType = 'combined';
        multScen.includeNomScen=true;
    case "impScen_permuted5"
        multScen = matRad_multScen(ct,'impScen');
        multScen.wcFactor=wcFactor;
        multScen.numOfShiftScen = [5 5 5];
        multScen.shiftSD = shiftSD;
        multScen.shiftGenType = 'equidistant';
        multScen.shiftCombType='permuted';
        multScen.numOfRangeShiftScen=124;
        multScen.rangeRelSD=0;
        multScen.rangeAbsSD=0;
        multScen.scenCombType = 'combined';
        multScen.includeNomScen=true;
    case "impScen_permuted7"
        multScen = matRad_multScen(ct,'impScen');
        multScen.wcFactor=wcFactor;
        multScen.numOfShiftScen = [7 7 7];
        multScen.shiftSD = shiftSD;
        multScen.shiftGenType = 'equidistant';
        multScen.shiftCombType='permuted';
        multScen.numOfRangeShiftScen=342;
        multScen.rangeRelSD=0;
        multScen.rangeAbsSD=0;
        multScen.scenCombType = 'combined';
        multScen.includeNomScen=true;
    case "impScen_permuted_truncated5"
        multScen = matRad_multScen(ct,'impScen');
        multScen.wcFactor=wcFactor;
        multScen.numOfShiftScen = [5 5 5];
        multScen.shiftSD = shiftSD;
        multScen.shiftGenType = 'equidistant';
        multScen.shiftCombType='permuted_truncated';
        multScen.numOfRangeShiftScen=32;
        multScen.rangeRelSD=0;
        multScen.rangeAbsSD=0;
        multScen.scenCombType = 'combined';
        multScen.includeNomScen=true; 
    case "impScen_permuted_truncated7"
        multScen = matRad_multScen(ct,'impScen');
        multScen.wcFactor=wcFactor;
        multScen.numOfShiftScen = [7 7 7];
        multScen.shiftSD = shiftSD;
        multScen.shiftGenType = 'equidistant';
        multScen.shiftCombType='permuted_truncated';
        multScen.numOfRangeShiftScen=122;
        multScen.rangeRelSD=0;
        multScen.rangeAbsSD=0;
        multScen.scenCombType = 'combined';
        multScen.includeNomScen=true;
    case "random"
        multScen = matRad_multScen(ct,'rndScen'); % 'impSamp' or 'wcSamp'
        multScen.probDist = 'equalProb';
        multScen.wcFactor=wcFactor;
        multScen.shiftSD = shiftSD;
        multScen.shiftGenType = 'sampled';
        multScen.shiftCombType = 'combined';
        multScen.numOfShiftScen = run_config.sampling_size  * ones(3,1);
        multScen.numOfRangeShiftScen = run_config.sampling_size;
        multScen.rangeRelSD=0;
        multScen.rangeAbsSD=0;
        multScen.scenCombType = 'combined';
        multScen.includeNomScen=true;
    case "random_truncated"
        multScen = matRad_multScen(ct,'rndScen'); % 'impSamp' or 'wcSamp'
        multScen.wcFactor=wcFactor;
        multScen.shiftSD = shiftSD;
        multScen.shiftGenType = 'sampled_truncated';
        multScen.shiftCombType = 'combined';
        multScen.numOfShiftScen = run_config.sampling_size * ones(3,1);
        multScen.numOfRangeShiftScen = run_config.sampling_size;
        multScen.rangeRelSD=0;
        multScen.rangeAbsSD=0;
        multScen.scenCombType = 'combined';
        multScen.includeNomScen=true;
    otherwise
        multScen = matRad_multScen(ct,'nomScen');
end

end

