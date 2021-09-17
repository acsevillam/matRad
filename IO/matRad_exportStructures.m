function matRad_exportStructures(cst,metadata)
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

if nargin<2
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
    if ~(isfolder('strExport'))
        mkdir(matRad_cfg.matRadRoot, 'strExport');
    end
    
    folderPath = [matRad_cfg.matRadRoot filesep 'strExport' filesep];
    
    [num_Struct, ~] = size(cst);
    
    %Create a file for each beam
    for i = 1:num_Struct
        
        if isempty(cst{i,4}) == false
            
            %Set a filename for i-th beam file
            filename_ith = [cst{i,2}];
            
            %Add column headers
            header_ith = header;
            header_ith = header_addComment(header_ith,'voxelID');
            
            data = cst{i,4};
            
            if strcmp(metadata.extension,'txt')
                
                %Write Header to file with the separating blank line to i-th beam
                fileHandle = fopen([folderPath filename_ith '.' metadata.extension],'w');
                fprintf(fileHandle,'%s\n',header_ith);
                
                %Append data to file to i-th beam
                writecell(data,[folderPath filename_ith '.' metadata.extension]);
                %dlmwrite(filename_ith+"."+metadata.extension,data{1,1},'delimiter',metadata.delimiter,'-append');
                
                fclose(fileHandle);
                
            elseif strcmp(metadata.extension,'bin')
                
                %Append data to file to i-th beam
                fileHandle = fopen([folderPath filename_ith '.' metadata.extension],'w');
                fwrite(fileHandle,uint32(data{1,1}),'uint32')
                fclose(fileHandle);
                
                %Write an additional header file
                headerHandle = fopen([folderPath filename_ith "_header.txt"],'w');
                fprintf(headerHandle,'%s\n',header_ith);
                fclose(headerHandle);
                
            end
            
        end
        
        fprintf(1,'structures exported successfully into %s.\n',strcat(folderPath,filename_ith,'.',metadata.extension));
        
    end
    
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