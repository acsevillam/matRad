function matRad_exportScenProb(filename,pln,metadata)
% matRad physical dose writer
%
% call
%   matRad_exportDij(filename,pln,...
%                    additionalFields,additionalKeyValuePairs)
%
% input
%   filename:   full output path, including the file extension
%   pln:        matRad pln struct
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

%% Setup scen_prob

scenProb_header = [];

%add matRad specific comment
scenProb_header = header_addComment(scenProb_header,'Created With matRad - An open source multi-modality radiation treatment planning sytem');

%% Write File
try
    
    %Set up parent export folder and full file path
    if ~(isfolder('dijExport'))
        mkdir(matRad_cfg.matRadRoot, 'dijExport');
    end
    
    folderPath = [matRad_cfg.matRadRoot filesep 'dijExport' filesep];
    
    scenProb = sprintf('#%s \t %s \t %s \t %s \t %s \n','scenNo','probability','dx','dy','dz');
    
    %Add info about scenario probabilities
    for i = 1:pln.multScen.totNumShiftScen
        scenProb = sprintf('%s%d \t %d \t %d \t %d \t %d \n',scenProb,i,pln.multScen.scenProb(i),pln.multScen.isoShift(i,1),pln.multScen.isoShift(i,2),pln.multScen.isoShift(i,3));
    end
    
    scenProbHandle = fopen([folderPath filename],'w');
    fprintf(scenProbHandle,'%s',scenProb_header);
    fprintf(scenProbHandle,'%s',scenProb);
    fclose(scenProbHandle);
    
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

fprintf(1,'Scenario probabilities exported successfully into %s.\n',strcat(folderPath,filename));

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