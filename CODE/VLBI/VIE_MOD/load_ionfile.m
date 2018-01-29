% ************************************************************************
%   Description:
%   Loads the external ionosperic delays. The delays have to be calculated
%   in after a previous run, using the azel file.
% 
%   References: 
%
%   Input:										
% 
%   Output:
% 
%   External calls: 	   
%       
%   Coded for VieVS: 
%   22 May 2012 by Lucia Plank
%
%   Revision: 
%   2016-10-17, A. Hellerschmied: Search path for .ion files set to: ionPath='../ION/FILES/';
%   
% ************************************************************************
function [iondata,ionFileFoundLog] = load_ionfile (parameter,session)


    ionPath='../ION/FILES/';

    ionSubfolder=parameter.vie_init.ionoFolder;

    % if the ext iono file is in main iono path
    if strcmp(ionSubfolder, ' ')
        ionSubfolder='/';
    end

    ionFolder=[ionPath, ionSubfolder, '/'];

    % define trpfile
    ionFile=[ionFolder, session, '.ion'];

    % if it does not exist, try to find it somewhere else (in
    % any subfolder)
    if ~exist(ionFile, 'file')
        ionFoundFiles=dirr([ionPath, '*', session(1:9)]);

        % if nothing found:
        if isempty(ionFoundFiles)
            fprintf('ion file not found\nNo ionospheric correction is applied\n\n');
            ionFileFoundLog=0;
        else
            % create new var just for search
            tempIonFoundFiles=ionFoundFiles(1); %(1)
            ionSearchSubFolder='/';

            while isstruct(tempIonFoundFiles)
                ionSearchSubFolder=[ionSearchSubFolder, tempIonFoundFiles.name, '/'];
                tempIonFoundFiles=tempIonFoundFiles.isdir;
            end
            % delete '/' at end
            ionSearchSubFolder(end)=[];

            % define new azel file
            ionFile=[ionPath, ionSearchSubFolder];
            fprintf('.ion file not found in specified folder.\nInstead following file found:\n''%s'' - this is being used\n', ionFile);
            ionFileFoundLog=1;
        end
    else
        fprintf('Reading external ion. file: %s\n', ionFile);
        ionFileFoundLog=1;
    end

    % if file found -> load and read
    if ionFileFoundLog==1
        % open
        fidIon=fopen(ionFile);

% frewind(fidIon)
        % read data
%         iondata = textscan(fidIon,'O %s %5f %4.0d.%2.0d.%2.0d-%2.0d:%2.0d:%4.1f %8s  %9f %8f  %6f %5f  %15f', 'CommentStyle', '#', 'delimiter', '||');
         iondata = textscan(fidIon,'O %s %d %4.0d.%2.0d.%2.0d-%2.0d:%2.0d:%4.1f %s %f %f  %f %f  %f', 'CommentStyle', '#');

        fclose(fidIon);
    end