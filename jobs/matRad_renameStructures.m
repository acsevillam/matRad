function [cst] = matRad_renameStructures(cst,run_config)

description=run_config.description;

switch description
    case 'prostate'

        for  it = size(cst,1):-1:1
            switch cst{it,2}
                case 'Skin'
                    cst{it,2}='BODY';
                case 'PTV'
                    cst{it,2}='PTV';
                case 'PROSTATA'
                    cst{it,2}='CTV';
                case 'RECTO'
                    cst{it,2}='RECTUM';
                case 'VEJIGA'
                    cst{it,2}='BLADDER';
                case 'BULBO'
                    cst{it,2}='BULB';
                case 'CFI'
                    cst{it,2}='LEFT REMORAL HEAD';
                case 'CFD'
                    cst{it,2}='RIGHT REMORAL HEAD';
                otherwise
                    fprintf(' %s is an empty structure. \n',cst{it,2});
                    fprintf('Deleting %s structure. \n',cst{it,2});
                    cst(it,:) = [];
            end
        end

    case 'breast'

        skin_flag = true;
        ixPTV=0;
        ixCTV=0;

        for  it = size(cst,1):-1:1
            switch cst{it,2}
                case 'Skin'
                    cst{it,2}='BODY';
                case 'Piel'
                    skin_flag = false;
                    cst{it,2}='SKIN';
                case {'PTV','PTV M'}
                    ixPTV=it;
                    cst{it,2}='PTV';
                case {'SENO IZQUIERDO','CTV'}
                    ixCTV=it;
                    cst{it,2}='CTV';
                case 'CORAZON'
                    cst{it,2}='HEART';
                case 'PULMON IZQUIERDO'
                    cst{it,2}='LEFT LUNG';
                case 'PULMON DERECHO'
                    cst{it,2}='RIGTH LUNG';
                case 'SENO CONTRALATERAL'
                    cst{it,2}='CONTRALATERAL BREAST';
                case 'MEDULA'
                    cst{it,2}='SPINAL CORD';
                otherwise
                    fprintf(' %s is an empty structure. \n',cst{it,2});
                    fprintf('Deleting %s structure. \n',cst{it,2});
                    cst(it,:) = [];
            end

            if(skin_flag && ixPTV~=0 && ixCTV~=0)
                metadata.name='SKIN';
                metadata.type='OAR';
                metadata.visibleColor=[1,0.501960784313726,1];
                [cst,~] = matRad_createSkin(ixPTV,ixCTV,cst,metadata);
                clear metadata;
            end

        end

end

for  it = 1:size(cst,1)
    cst{it,1}=it-1;
end

end

