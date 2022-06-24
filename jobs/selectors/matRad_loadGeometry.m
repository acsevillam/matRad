function [ct,cst] = matRad_loadGeometry(run_config)

    patient=run_config.description;
    resolution=run_config.resolution;
    
    switch patient
        case 'prostate'
            if(resolution=="3x3x3")
                load('patient1_3x3x3mm.mat');
            end

            if(resolution=="5x5x5")
                load('patient1_5x5x5mm.mat');
            end  
        case 'breast'
            if(resolution=="3x3x3")
                load('patient2_3x3x3mm.mat');
            end

            if(resolution=="5x5x5")
                load('patient_2_multiScen_contourned.mat');
            end 
        case 'H&N'
            if(resolution=="5x5x5")
                load('patient4_5x5x5mm.mat');
            end 
    end

end

