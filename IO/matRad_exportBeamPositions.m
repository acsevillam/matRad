function matRad_exportBeamPositions(pln,stf,metadata)
% matRad structures writer
%
% call
%   matRad_exportStructures(cst,...
%                    additionalFields,additionalKeyValuePairs)
%
% input
%   cst:        matRad cst struct
%   metadata:   struct of metadata
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

if nargin<3
    metadata = struct();
end

%% Prepare Metadata

if ~isfield(metadata,'delimiter')
    metadata.delimiter = '\t'; %Default delimiter
end

if ~isfield(metadata,'extension')
    metadata.extension = 'txt'; %Default fileType
end

%% Setup Header

header = sprintf('# %s %s\n',metadata.extension,'file');

%add matRad specific comment
header = header_addComment(header,'Created With matRad - An open source multi-modality radiation treatment planning sytem');

%% Write File
try
    
    %Set up parent export folder and full file path
    if ~(isfolder('beamExport'))
        mkdir(matRad_cfg.matRadRoot, 'beamExport');
    end
    
    folderPath = [matRad_cfg.matRadRoot filesep 'beamExport' filesep];
    
    filename = 'beam_pos.txt';
    fileHandle1 = fopen([folderPath filename],'w');
    fileHandle2 = fopen([folderPath filename '_header.txt'],'w');
    
    for i = 1:pln.propStf.numOfBeams
        fprintf(fileHandle2,'%i\n',stf(i).numOfRays);
        for j = 1:stf(i).numOfRays
            fprintf(fileHandle1,'%i\t%i\t%i\n',stf(i).ray(j).rayPos_bev(1),stf(i).ray(j).rayPos_bev(2),stf(i).ray(j).rayPos_bev(3));
        end
    end
    
    fprintf(1,'beam positions exported successfully into %s.\n',strcat(folderPath,filename,'.',metadata.extension));
    
catch MExc
    %if something failed while writing, close all files and display error
    fclose('all');
    fprintf(2,'File %s could not be written!\n',filename);
    if(matRad_cfg.isOctave)
        error(MExc);
    else
        throw(MExc);
    end
end

%Used to add comments to the header
    function newHeader = header_addComment(header,comment)
        newHeader = sprintf('%s# %s\n',header,comment);
    end

%Used to add int fields to the header
    function newHeader = header_addIntField(header,fieldName,fieldValue)
        newHeader = sprintf('%s# %s: %d\n',header,fieldName,fieldValue);
    end

%Used to add string fields to the header
    function newHeader = header_addStringField(header,fieldName,fieldValue)
        newHeader = sprintf('%s# %s: %s\n',header,fieldName,fieldValue);
    end

end

