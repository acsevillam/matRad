function matRad_exportStructures(filename,cst,metadata)
% matRad structures writer
%
% call
%   matRad_exportStructures(filename,cst,...
%                    additionalFields,additionalKeyValuePairs)
%
% input
%   filename:   full output path, including the file extension
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
if nargin<4
    metadata = struct();
end


%% Prepare Metadata

if ~isfield(metadata,'delimiter')
    metadata.delimiter = '\t'; %Default delimiter
end

if ~isfield(metadata,'numScen')
    metadata.numScen = 1; %Default scenario
end

if ~isfield(metadata,'extension')
    
    lastdot_pos = find(filename == '.', 1, 'last');
    extension = filename(lastdot_pos+1:end);
    
    if strcmp(extension,'txt') || strcmp(extension,'bin')
        metadata.extension = extension; %Default fileType
    else
        metadata.extension = 'txt'; %Default fileType
    end
    
end

%% Setup Header

header = sprintf('# %s %s\n',metadata.extension,'file');

%add matRad specific comment
header = header_addComment(header,'Created With matRad - An open source multi-modality radiation treatment planning sytem');

%% Write File
try
        
    [num_Struct, ~] = size(cst);

    %Create a file for each beam
    for i = 1:num_Struct

        if isempty(cst{i,4}) == false

            %Set a filename for i-th beam file
            lastdot_pos = find(filename == '.', 1, 'last');

            filename_ith = filename(1:lastdot_pos-1);
            filename_ith = filename_ith+"_"+cst{i,2};

            %Add column headers
            header_ith = header;
            header_ith = header_addComment(header_ith,'voxelID');

            data = cst{i,4};
            
            if strcmp(metadata.extension,'txt')

                %Write Header to file with the separating blank line to i-th beam
                fileHandle = fopen(filename_ith+"."+metadata.extension,'w');
                fprintf(fileHandle,'%s\n',header_ith);

                %Append data to file to i-th beam
                writecell(data,filename_ith+"."+metadata.extension);
                %dlmwrite(filename_ith+"."+metadata.extension,data{1,1},'delimiter',metadata.delimiter,'-append');

                fclose(fileHandle);

            elseif strcmp(metadata.extension,'bin')

                %Append data to file to i-th beam
                fileHandle = fopen(filename_ith+"."+metadata.extension,'w');
                fwrite(fileHandle,uint32(data{1,1}),'uint32')
                fclose(fileHandle);

                %Write an additional header file
                headerHandle = fopen(filename_ith+"_header.txt",'w');
                fprintf(headerHandle,'%s\n',header_ith);
                fclose(headerHandle);

            end

        end

    end
    
catch MExc
    fclose('all');
    error(sprintf('File %s could not be written!\n%s',filename,getReport(MExc)));
    
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