% ************************************************************************
%   Description:
%   Loads the external tropospheric delays. The delays have to be calculated
%   in a previous run, using the azel file.
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
% 
%   07 Apr 2014 by A. Hofmeister: include checking of existing trp-files which don't include the
%                                    NGS-file version (_N004) in its session name (e.g. when using
%                                    GSFC files) in case that normal trp-files with version included
%                                    in session name are not found,
%                                    report if trp-files without version identifier will be used
%   28 Jul 2014 by A. Hofmeister: add query about trp-format of input file (old (Matthias)
%                                    format or new format),
%                                    correct comments,
%                                    add comments to structure the function sections
%   29 Jul 2014 by A. Hofmeister: programming of read in of new format
%   30 Jul 2014 by A. Hofmeister: add comments
%   16 Dec 2014 by D. Landskron:  partial Gn and Ge deleted (not used)
%   22 Dec 2015 by A. Hofmeister: remove old file format specifications for reading the trp-file,
%                                    remove check of header in order that trp-files of TU Wien
%                                    and GSFC can be read in with the same format specifications,
%                                    make data line recognition more robust,
%                                    correct comments,
%                                    remove unnecessary old output parameter "vahabsRaytracingFiles"
%   06 Oct 2016 by A. Hofmeister: set trpdata=[] if no trp file is found since trpdata is a
%                                    needed output variable
%                                    remove message about backup for trpdata, which is displayed
%                                    later by vie_mod.m
%   20 Jan 2017 by D. Landskron: enabled also for ray-tracing .trp-files
%   10 Mar 2017 by D. Landskron: raytr data also enabled from yearly subdivided folders
%   13 Sep 2017 by D. Landskron: 'tropSource' shifted into 'vie_init' 
%   09 Jul 2018 by D. Landskron: also enabled for input of .radiate files
%   01 Aug 2018 by D. Landskron: .radiate files made priority input
%
%
function [trpdata,trpFileFoundLog] = load_trpfile (parameter,session)


%% Get .trp-file or .radiate-file for the desired session
% This section is used to find the correct .trp-files or .radiate-files for the desired session.
% In case that the predefined (user-defined) path to file does not contain the file, it will be
% tried to exclude the ngs-file version (e.g. "_N004"). If the file is still not found all
% subfolders (and sub-...-subfolders) in the main trp-path will be searched to find a trp-file
% matching with the session name (except for the ngs-file version).


if ~strcmp(parameter.vie_init.tropSource.name,'raytr')
    error('Something is wrong with the usage of raytracing files!')
end
    

% define the raytracing file
trpFolder = '../TRP/RAYTRACING_DATA/';
if exist([trpFolder,parameter.year,'/'])==7
    trpFile=[trpFolder,parameter.year,'/' session, '.trp'];
    radiateFile=[trpFolder,parameter.year,'/' session, '.radiate'];
else
    trpFile=[trpFolder, session, '.trp'];
    radiateFile=[trpFolder, session, '.radiate'];
end

   
% if both .radiate-file and .trp file of a session are available, then the .radiate-file is read, because it's more accurate
if exist(radiateFile, 'file')
    radiateFileFoundLog=1;
    trpFileFoundLog=0;
elseif exist(trpFile, 'file')
    radiateFileFoundLog=0;
    trpFileFoundLog=1;
elseif ~exist(trpFile, 'file') 
    error('No ray-tracing file (.trp) available for this session! Specify another source for the tropospheric delays.')
end


% This section reads in the trp-file in case one has been found.
% A new variable "trpdata" will be created, which must have the following structure for a correct
% further usage in VieVS:
%
% Content of variable "trpdata":
%   trpdata is a cell array of dimension 1x15
%   content of the cells
%       cell 1:  source name
%       cell 2:  scan number, info not used in VieVS
%       cell 3:  year
%       cell 4:  month
%       cell 5:  day
%       cell 6:  hour
%       cell 7:  minute
%       cell 8:  seconds
%       cell 9:  site ID
%       cell 10: azimuth in [°], info not used in VieVS
%       cell 11: outgoing elevation angle in [°], info not used in VieVS
%       cell 12: pressure in [hPa]
%       cell 13: temperature in [°C]
%       cell 14: slant total delay including geometric bending effect in [sec]
%       cell 15: wet mapping function
%
% Note: The additional data columns contained in the .trp-file and .radiate
% file are not read in



%% (a) Read in the .radiate-file
% if both .trp-file and .radiate file of a session are available, then the .radiate-file is read, because it's more accurate


% if file found -> load and read
if radiateFileFoundLog==1
    
    c = 299792458;   % speed of light in [m/s]
    
    % open the file
    fidRadiate=fopen(radiateFile);
    
    % get the current line
    curr_line=fgetl(fidRadiate);
    curr_line=fgetl(fidRadiate);
    if ~strcmpi(curr_line,'% RADIATE format v 2.0')
        error('Something is wrong with the .radiate file!');
    end
    fclose(fidRadiate);
    
    % open the file again and read all data
    fidRadiate=fopen(radiateFile);    
    radiate_data = textscan(fidRadiate,'%f%f%f%f%f%f%f%s%f%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','CommentStyle','%');
    
    mjd = radiate_data{2};
    [~, month, day, ~, ~, ~] = mjd2date(mjd);
    
    % define the correct columns
    trpdata{1,1} = radiate_data{11};           % source name
    trpdata{1,2} = radiate_data{1};            % scan number
    trpdata{1,3} = radiate_data{3};            % year
    trpdata{1,4} = month;                      % month
    trpdata{1,5} = day;                        % day
    trpdata{1,6} = radiate_data{5};            % hour
    trpdata{1,7} = radiate_data{6};            % min
    trpdata{1,8} = radiate_data{7};            % sec
    trpdata{1,9} = radiate_data{8};            % station
    trpdata{1,10} = radiate_data{9}*180/pi;    % azimuth angle in [°]
    trpdata{1,11} = radiate_data{10}*180/pi;   % outgoing elevation angle in [°]
    trpdata{1,12} = radiate_data{13};          % pressure in [hPa]
    trpdata{1,13} = radiate_data{12};          % temperature in [°C]
    trpdata{1,14} = radiate_data{18}/c;        % slant total delay including geometric bending effect in [sec]  
    trpdata{1,15} = radiate_data{26};          % wet mapping function
    
    % add blanks to the source names and station names
    trpdata{1} = pad(trpdata{1},8);
    trpdata{9} = pad(trpdata{9},8);
    
    % close the file
    fclose(fidRadiate);
     
end



%% (b) Read in the trp-file
% this is only read in if no .radiate-file is available


% if file found -> load and read
if trpFileFoundLog==1
    
    % open the file
    fidTrp=fopen(trpFile);
    
    % set variable for determining the index of the current O-record 
    ind=0;
    
    % loop till end of file
    while ~feof(fidTrp)
        
        % get the current line
        curr_line=fgetl(fidTrp);
        
        % check if current line contains an O-record, which contains the necessary trp-data
        % --> read in of data
        % note: Add test if lenght of current line is larger than two characters as it will be
        %       needed for the test if the current line is a data line.
        %       Test for 'O  ' including two blanks following the 'O' as specified in the trp-format
        %       to avoid that a modified header at the beginning of the file or at the end of the
        %       file is accidently interpreted as data line if it starts with 'O'.
        if length(curr_line) > 2
            if strcmp(curr_line(1:3), 'O  ')
            
                % raise index of current record
                ind=ind+1;

                % extract the different values from the record line (observation and ray-tracing
                % data)
                % use the trp-format specification (see trp-file) to find the necessary span of
                % characters describing each value

                trpdata{1,1}{ind,1} = curr_line(13:20);                  % source name
                trpdata{1,2}(ind,1) = str2double(curr_line(4:8));        % scan number
                trpdata{1,3}(ind,1) = str2double(curr_line(26:29));      % year
                trpdata{1,4}(ind,1) = str2double(curr_line(31:32));      % month
                trpdata{1,5}(ind,1) = str2double(curr_line(34:35));      % day
                trpdata{1,6}(ind,1) = str2double(curr_line(37:38));      % hour
                trpdata{1,7}(ind,1) = str2double(curr_line(40:41));      % minute
                trpdata{1,8}(ind,1) = str2double(curr_line(43:46));      % seconds
                trpdata{1,9}{ind,1} = curr_line(49:56);                  % station
                trpdata{1,10}(ind,1) = str2double(curr_line(59:67));     % azimuth angle in [°]
                trpdata{1,11}(ind,1) = str2double(curr_line(69:76));     % outgoing elevation angle in [°]
                trpdata{1,12}(ind,1) = str2double(curr_line(79:84));     % pressure in [hPa]
                trpdata{1,13}(ind,1) = str2double(curr_line(86:90));     % temperature in [°C]
                trpdata{1,14}(ind,1) = str2double(curr_line(93:107));    % slant total delay including geometric bending effect in [sec]
                trpdata{1,15}(ind,1) = str2double(curr_line(109:123));   % wet mapping function
                
            end
        end
        
    end % end of loop over lines in file
    
    % close the file
    fclose(fidTrp);
    
end




end