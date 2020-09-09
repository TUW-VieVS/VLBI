% This function creates external ionospheric files. The model used (GNSS,
% GNSS+Altimetry, GNSS+Altimetry+F/C, CODE, IGS) can be chosen via GUI
% start_createExtIonoFilesGUI. This is also the main GUI to do this task
% (creating ionospheric files).
%
% INPUT
%       sessionName     The name of the VLBI session for which the iono
%                       external file should be created.
%       ionSubDir       The subfolder of '\VieVS\ION\?\' where the data
%                       should be stored#
%       azelFile        Full path and file of the needed azel file
%                       (containing azimuth and elevation information).
%       ionoModel       Defines which model to be used (one of:
%                       - 'CODE'
%                       - 'IGS'
%
% OUTPUT (optional)
%       errorMessage    Integer specifing an error caught:
%                       0 | no error
%                       1 | name of subfolder is not allowed
%                       2 | azelFile could not be created (probably need to
%                           process session in VieVS)
%                       3 | IONO map does not exist -> was downloaded (need
%                           to be unzipped)
%                       4 | GNSS/GNSS+Altimetry/GNSS+Altimetry+F/C not yet
%                           implemented
%                       5 | antenna struct not found (there was search
%                           performed
%                       6 | station (in azel) was not found in antenna struct
%                       7 |
%
% CHANGES
%   2016-10-14, A. Hellerschmied: - Reviewed and adjusted function for VieVS 3.0
%                                 - "totDelStr(end-2) = []" ("-" from exponential representation of the total ion. delay NOT removed any more) commented
%   2017-04-04, A. Hellerschmied: Ref. frequ. can now be changed in CW
%   2017-08-31, A. Hellerschmied: Bug fix
%
%**************************************************************************


function varargout=createExtIonoFiles(sessionName, ionSubDir, azelFile, ionoModel)

% ##### Options: #####
ref_freq_Hz = 8.4e9; % X-band (default)
% ref_freq_Hz = 1567e6; % L1

% Status Msg.:
fprintf('Creating ion. file for session: %s\n', sessionName);
fprintf('Reference frequency =  %5.2f MHz (default for X-band)\n', ref_freq_Hz * 1e-6);


% Change ref. frequ. 
flag_input_accepted = false;
while(~flag_input_accepted)
    input_str = input(' => Please enter "y" to take the default ref. frequ., or enter an alternative frequ. [MHz]:', 's');
    switch(input_str)
        case 'y'
            fprintf('Reference frequency =  %5.2f MHz\n', ref_freq_Hz * 1e-6);
            flag_input_accepted = true;
        otherwise
            [ref_freq_MHz, status] = str2num(input_str);
            if status ~= 1
                fprintf('       - ERROR: Invalid input: Only numbers allowed (or "y")!\n')
            else
                ref_freq_Hz = ref_freq_MHz * 1e6;
                fprintf('Reference frequency =  %5.2f MHz\n', ref_freq_MHz);
                flag_input_accepted = true;
            end
    end
end % while(~flag_input_accepted)


%~~~~~~~~~~~~~~~~~~~~~~~~ DEFINITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ##### Define filepaths #####
% Output Ion. file:
fidOutPath      = ['../ION/FILES/', ionSubDir, '/'];
fidOutFile      = [sessionName, '.ion'];
% Ion. maps:
ionoPath        = repmat('../ION/MAPS/', 2, 1);
% Antenna structure in LEVEL1 directory (output from VIE_MOD)
antPath         = '../DATA/LEVEL1/';
antennaStruct   = [antPath, sessionName, '_antenna.mat'];

% ##### add path to auxiliary code #####
addpath('../ION/PROGRAM/auxiliary/');

% ##### Get constants #####
constants
global c


% ##### Init.: #####
level1_subdir_str = '';


%~~~~~~~~~~~~~~~~~~~~~~~~~~ FOLDERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% check name of subfolders
forbiddenSubfolders = {'PROGRAM', 'MAPS'};
if sum(ismember(forbiddenSubfolders, ionSubDir)) ~= 0
    errorMessage = 1;
    varargout(1) = {errorMessage};
    return;
end

% create subfolder if not exist
if ~exist(fidOutPath, 'dir')
    mkdir(fidOutPath);
end

% folder of ionosphere map
if ~exist(ionoPath(1,:), 'dir')
    mkdir(ionoPath(1,:));
end

% card/map folder ('..\CODE\')
ionoPath = repmat([ionoPath(1,:), ionoModel, '/'], 2, 1);
if ~exist(ionoPath(1,:), 'dir')
    mkdir(ionoPath(1,:));
end


%~~~~~~~~~~~~~~~~~~~~~~~ LOAD ANTENNA STRUCT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fprintf(' => Load antenna structure from LEVEL1 dir.: %s\n', [sessionName, '_antenna.mat']);

if ~exist(antennaStruct, 'file')
    
    % #### Get directories, where the required antenna file is located: ####
    antFoundFiles = dirr([antPath, '*', sessionName, '_antenna.mat']);
    
    % Check, if the antenna file was found:
    if isempty(antFoundFiles)
        errorMessage = 5;
        varargout(1) = {errorMessage};
        return;
    end
    
    % #### In case of multiple results => User has to select the correct directory: ####
    % Only one dir.:
    if length(antFoundFiles) == 1
        fprintf('   - antenna struct. found in dir.: %s\n', [antPath, antFoundFiles.name])
        tempAntFoundFiles = antFoundFiles(1);
        level1_subdir_str = tempAntFoundFiles.name;
        
        % Multiple directories:
    else
        fprintf('   - antenna struct. found in multiple directories.: %s\n', [antPath, antFoundFiles.name])
        for i_tmp = 1 : length(antFoundFiles)
            fprintf('      - %2d : %s\n',i_tmp, [antPath, antFoundFiles(i_tmp).name])
        end
        
        % User has to select the LEVEL1 sub.-dir.:
        flag_input_accepted = false;
        while(~flag_input_accepted)
            input_str = input('       => Please select the correct LEVEL1 sub-dir. (enter number):', 's');
            [selected_dir, status] = str2num(input_str);
            % ### Check input: ###
            % Is a number?
            if status ~= 1
                fprintf('       - ERROR: Invalid input: Only numbers allowed!\n')
            else
                if (selected_dir < 1) || (selected_dir > length(antFoundFiles))
                    fprintf('       - ERROR: Invalid input: Enter a number between 1 and %d!\n', length(antFoundFiles));
                else
                    flag_input_accepted = true;
                    fprintf('       - Selected sub.-dir.: %s\n', antFoundFiles(selected_dir).name)
                end
            end
        end % while(~flag_input_accepted)
        tempAntFoundFiles = antFoundFiles(selected_dir);
        level1_subdir_str = tempAntFoundFiles.name;
    end
    
    % create new var just for search
    antSearchSubFolder='/';
    
    while isstruct(tempAntFoundFiles)
        antSearchSubFolder = [antSearchSubFolder, tempAntFoundFiles.name, '/'];
        tempAntFoundFiles = tempAntFoundFiles.isdir;
    end
    
    % delete '\' at end
    antSearchSubFolder(end)=[];
    
    % define new antenna file
    antennaStruct=[antPath, antSearchSubFolder];
    
else
    antSearchSubFolder='';
end
% load antenna struct
load(antennaStruct);



%~~~~~~~~~~~~~~~~~~~~~~~ READ AZEL FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fprintf(' => Load AZEL file: %s\n', azelFile);


% if it does not exist -> try to create it
if ~exist(azelFile, 'file')
    fprintf('   - AZEL file does not exist => Try to create it!\n')
    
    % Get path:
    str_ind = strfind(azelFile, '/');
    path_azel_file_str = azelFile(1 : str_ind(end-1));
    
    % Get year:
    year_str = azelFile(str_ind(end-1)+1 : str_ind(end)-1);
    
    % create process_list
    process_list = [year_str, '/', sessionName];
    
    % LEVEL1 sub.-dir.:
    subdir_azel_str = level1_subdir_str;
    
    % call azel_out.m to create azel file:
    azel_out(process_list, subdir_azel_str, path_azel_file_str)
    
end

% if file now exists (after creating)
if exist(azelFile, 'file')
    
    % Read AZEL file:
    fidAzel     = fopen(azelFile);
    azel_data   = textscan(fidAzel, '%4.0f %11.5f %4.0f %3.0f %2.0f %2.0f %5.2f %8s %20.15f %20.15f %8s %6.2f %7.2f %*[^\n]', 'CommentStyle', '%');
    fclose(fidAzel);
    
    % get number of lines of azel file:
    nObs = size(azel_data{1},1);
    
    % make all statNames 8 digits long
    for k = 1 : nObs
        
        curLength = length(azel_data{8}{k}); % length of current station
        if curLength < 8
            azel_data{8}{k}(curLength+1 : 8) = ' ';
        end
    end
else % if azel file was not found though creating
    errorMessage=2;
    varargout(1)={errorMessage};
    return;
end


%~~~~~~~~~~~~~~~~~~~~~~~ DATE CONVERSIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Year (for both days)
curYr           = azel_data{3}(1);
curYrStr        = num2str(curYr);
curYr(2,1)      = azel_data{3}(end);
curYrStr(2,:)   = num2str(curYr(2,1));

% Day of year (for both days)
yyDoySecod = [azel_data{3}(1), azel_data{4}(1), 0; azel_data{3}(end), azel_data{4}(end), 0];

% mjd (for both days)
curMjd      = floor(azel_data{2}(1));
curMjd(2,1) = floor(azel_data{2}(end));

% Date (for both days)
[~, month1st day1st] = mjd2date(curMjd(1,1));
[~, month2nd day2nd] = mjd2date(curMjd(2,1));



%~~~~~~~~~~~~~~~~~~~~~~~~ DOWNLOAD IONEX FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% define variable for not found ionex
% currentfolder = pwd; %gets currently working directory

for k = 1 : 2 % loop over 2 days
    if strcmp(ionoModel, 'CODE') == 1 % if CODE map is being used
        
        curIonoPath = [ionoPath(k,:), curYrStr(k,:), '/'];
        if ~exist(curIonoPath, 'dir') % if directory does not exist, create it
            mkdir(curIonoPath);
        end
        ionoFilename{k,1} = ['CODG', sprintf('%03.0f', yyDoySecod(k,2)), '0.', curYrStr(k,3:4), 'I'];
        ionoFile{k,1} = [curIonoPath, ionoFilename{k,1}];
        
        if ~exist(ionoFile{k,1}, 'file') % if the ionex file is not existing unzipped
            fprintf('CODE TEC map file does not exist:            (%s)\n', ionoFile{k,1});
            zippedFile = [curIonoPath, ionoFilename{k,1}, '.Z'];
            if ~exist(zippedFile, 'file') % if zipped file does not exist in folder -> download
                fprintf('Compressed IGS TEC map file does not exist: (%s)\n', zippedFile);
%                 url = ['ftp://ftp.unibe.ch/aiub/CODE/', num2str(curYr(k)), '/', ionoFilename{k,1}, '.Z'];
                url = ['http://ftp.aiub.unibe.ch/CODE/', num2str(curYr(k)), '/', ionoFilename{k,1}, '.Z'];

                urlwrite(url, zippedFile);
                fprintf(' ...finished downloading.\n');
            end
            
            % write message to user to unzip manually (and return afterwards)
            fprintf('Zipped file (''%s'') was downloaded or already exists in\n''%s''\nPlease extract manually (to same folder)!\n\n', [ionoFilename{k,1}, '.Z'], [ionoPath(k,:), num2str(curYr(k)), '/'] );
            
             % User has to select the LEVEL1 sub.-dir.:
            flag_input_accepted = false;
            while(~flag_input_accepted)
                input_str = input(' => Please enter "y" after uncrompressing the file to continue:', 's');
                % ### Check input: ###
                switch(input_str)
                    case 'y'
                        flag_input_accepted = true;
                    otherwise
                        fprintf('    - Error: Invalid input!')
                end
            end % while(~flag_input_accepted)
%             cd(curIonoPath)
%             unix(strcat(ionoFilename{k,1},'.Z uncompress ')); %uncompress
%             delete(strcat(ionoFilename{k,1},'.Z'));
%             cd(currentfolder)
            
        end
    else % if IGS map is being used
        
        curIonoPath=[ionoPath(k,:), curYrStr(k,:), '/'];
        if ~exist(curIonoPath, 'dir')  % if directory does not exist, create it
            mkdir(curIonoPath);
        end
        ionoFilename{k,1}=['igsg', sprintf('%03.0f', yyDoySecod(k,2)), '0.', curYrStr(k,3:4), 'i'];
        ionoFile{k,1}=[curIonoPath, ionoFilename{k,1}];
        
        if ~exist(ionoFile{k,1}, 'file') % if the ionex file is not existing unzipped
            fprintf('IGS TEC map file does not exist:            (%s)\n', ionoFile{k,1});
            zippedFile = [curIonoPath, ionoFilename{k,1}, '.Z'];
            if ~exist(zippedFile, 'file')  % if zipped file does not exist in folder -> download
                fprintf('Compressed IGS TEC map file does not exist: (%s)\n', zippedFile);
                url=['ftp://cddis.gsfc.nasa.gov/gps/products/ionex/', num2str(curYr(k)), '/', sprintf('%03.0f', yyDoySecod(k,2)), '/', ionoFilename{k,1}, '.Z'];
                fprintf(' => File will be downloaded from URL: %s\n', url);
                urlwrite(url, zippedFile);
                fprintf(' ...finished downloading.\n');
            end
            
            % write message to user to unzip manually (and return afterwards)
            fprintf('Zipped file (''%s'') was downloaded or already exists in\n''%s''\nPlease extract manually (to same folder)!\n\n', [ionoFilename{k,1}, '.Z'], [ionoPath(k,:), num2str(curYr(k)), '/'] );
            
            % User has to select the LEVEL1 sub.-dir.:
            flag_input_accepted = false;
            while(~flag_input_accepted)
                input_str = input(' => Please enter "y" after uncrompressing the file to continue:', 's');
                % ### Check input: ###
                switch(input_str)
                    case 'y'
                        flag_input_accepted = true;
                    otherwise
                        fprintf('    - Error: Invalid input!\n')
                end
            end % while(~flag_input_accepted)
%             cd(curIonoPath)
%             unix(strcat(ionoFilename{k,1},'.Z uncompress '));
%             delete(strcat(ionoFilename{k,1},'.Z'));
%             cd(currentfolder)
        end
    end
end



%~~~~~~~~~~~~~~~~~~~~~~~~~ READ IONEX FILES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flag_ionex_files_equal = strcmp(ionoFile{1,1}, ionoFile{2,1});

if flag_ionex_files_equal
    fprintf('Reading IONEX file: %s\n', ionoFile{1,1});
    io1 = ionex_read(ionoFile{1,1});
    io2 = io1;
else
    fprintf('Reading IONEX files: %s and %s\n', ionoFile{1,1}, ionoFile{2,1});
    io1 = ionex_read(ionoFile{1,1});
    io2 = ionex_read(ionoFile{2,1});
end


%~~~~~~~~~~~~~~~~~~~~ ADD COLUMNS TO azel_data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% new column number to write antenna index
newCol4azel = size(azel_data,2)+1;

% all unique stations
uniqueStat = unique(azel_data{8});

for k = 1 : size(uniqueStat, 1)
    % find this station in antenna struct, in vmf1Data cell
    curAntIdx = find(~cellfun(@isempty, strfind({antenna.name}, uniqueStat{k})));
    
    % if station was not found in antenna struct
    if isempty(curAntIdx)
        errorMessage = 6;
        varargout(1) = {errorMessage};
        return;
    end
    
    % get rows where antenna-struct-index should be written to (rows of current stations)
    rows2write = ~cellfun(@isempty, strfind(azel_data{8}, uniqueStat{k}));
    
    % write this one index to many lines (where thsi one station is)
    azel_data{newCol4azel}(rows2write,1) = curAntIdx;
    
    % write x,y,z to following columns of azel
    azel_data{newCol4azel+1}(rows2write,1) = antenna(curAntIdx).x;
    azel_data{newCol4azel+2}(rows2write,1) = antenna(curAntIdx).y;
    azel_data{newCol4azel+3}(rows2write,1) = antenna(curAntIdx).z;
    
    % write ellipsoidal lat, lon, h to azel
    [lat,lon,h] = xyz2ell([antenna(curAntIdx).x, antenna(curAntIdx).y, antenna(curAntIdx).z]);
    
    azel_data{newCol4azel+4}(rows2write,1) = lat;
    azel_data{newCol4azel+5}(rows2write,1) = lon;
    azel_data{newCol4azel+6}(rows2write,1) = h;
end


%~~~~~~~~~~~~~~~~~~~ CALCULATION OF DELAY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fprintf('Calculating delays.\n');

zenDist = pi/2-azel_data{10};

% calculate coefficients
b_tec_inter = -2*io1.radius*1000*cos(pi-zenDist);
c_tec_inter = -((io1.hgt1*1000)^2+2*(io1.radius*1000)*(io1.hgt1*1000));
d1          = (-b_tec_inter+sqrt(b_tec_inter.^2-4*c_tec_inter))/2;

dren(:,1) = d1.*cos(zenDist);
dren(:,2) = d1.*sin(zenDist).*sin(azel_data{9});
dren(:,3) = d1.*sin(zenDist).*cos(azel_data{9});

Incremento_xyz = ren2xyz(dren, azel_data{18}, azel_data{19});

xipp = azel_data{15} + Incremento_xyz(:,1);
yipp = azel_data{16} + Incremento_xyz(:,2);
zipp = azel_data{17} + Incremento_xyz(:,3);

% calculate latitude
latpp           = zeros(size(zipp,1),1); % preallocate
latpp(zipp>0)   = (pi/2-atan(sqrt(xipp(zipp>0).^2+yipp(zipp>0).^2)./zipp(zipp>0)));
latpp(zipp<=0)  = pi/2-(pi+atan(sqrt(xipp(zipp<=0).^2+yipp(zipp<=0).^2)./zipp(zipp<=0)));

% calculate longitude
lonpp                   = zeros(size(zipp,1),1); % preallocate
lonpp(xipp>0 & yipp>0)  = atan(yipp(xipp>0 & yipp>0)./xipp(xipp>0 & yipp>0));
lonpp(xipp>0 & yipp<0)  = 2*pi+atan(yipp(xipp>0 & yipp<0)./xipp(xipp>0 & yipp<0));
lonpp(xipp<=0)          = pi+atan(yipp(xipp<=0)./xipp(xipp<=0));

% Convert in Grad
latpp = latpp*180/pi;
lonpp = lonpp*180/pi;

lonpp(lonpp < -180) = lonpp(lonpp < -180)+360;
lonpp(lonpp > 180)  = lonpp(lonpp > 180)-360;
latpp(latpp < -90)  = -(180+latpp(latpp < -90));
latpp(latpp > 90)   = 180-latpp(latpp > 90);

% calculation of VTEC
VTEC = zeros(size(azel_data{1},1),1); % preallocate

% for all lines in azel (loop needed for interpolation function)
for k = 1 : size(azel_data{1}, 1)         
    if azel_data{4}(k) == azel_data{4}(1)
        VTEC(k,1) = interpolation(io1, latpp(k), lonpp(k), azel_data{5}(k), azel_data{6}(k), azel_data{7}(k)); %doy
    else
        VTEC(k,1) = interpolation(io2, latpp(k), lonpp(k), azel_data{5}(k), azel_data{6}(k), azel_data{7}(k));
    end
end

% calculation of M(zd), STEC
Mzd                                 = zeros(size(azel_data{1},1),1);
Mzd(azel_data{4}==azel_data{4}(1))  = 1./(sqrt(1-((io1.radius/(io1.radius+io1.hgt1)).*sin(zenDist(azel_data{4}==azel_data{4}(1)))).^2));
Mzd(azel_data{4}~=azel_data{4}(1))  = 1./(sqrt(1-((io2.radius/(io2.radius+io2.hgt1)).*sin(zenDist(azel_data{4}~=azel_data{4}(1)))).^2));
STEC                                = VTEC.*Mzd;

% calculate delay
ion =(40.28e16/ref_freq_Hz^2)*STEC/c; % [sec]


%~~~~~~~~~~~~~~~~~~~~~~~~ WRITE .ION FILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fprintf('Write ion. file: %s\n', [fidOutPath, fidOutFile]);

% get dates for observations
years2write     = zeros(size(azel_data{1},1),1);
months2write    = years2write;
days2write      = years2write;

years2write(azel_data{4}==azel_data{4}(1))  = curYr(1,1);
years2write(azel_data{4}~=azel_data{4}(1))  = curYr(2,1);
months2write(azel_data{4}==azel_data{4}(1)) = month1st;
months2write(azel_data{4}~=azel_data{4}(1)) = month2nd;
days2write(azel_data{4}==azel_data{4}(1))   = day1st;
days2write(azel_data{4}~=azel_data{4}(1))   = day2nd;

% open file for writing
fidOut=fopen([fidOutPath, fidOutFile], 'w');

% write header information
curClock = clock;
fprintf(fidOut, '# IONOSPHERIC DELAY\n#\n# This file was created on: %02.0f.%02.0f.%0.0f (%02.0f:%02.0f:%02.0f)\n#\n# Session: %s\n# IONEX-Map used: %s\n# Ref. frequ = %5.2f MHz\n#\n#\n', curClock(3), curClock(2), curClock(1), curClock(4), curClock(5), curClock(6), sessionName, ionoModel, ref_freq_Hz*1e-6);
fprintf(fidOut, '%s\n','#               X[m]		   Y[m]			  Z[m]	   geocLat(°) long(°)   ellHeight[m]');

% write stations
for stat = 1 : size(antenna, 2)
    [lat,lon,h] = xyz2ell([antenna(stat).x, antenna(stat).y, antenna(stat).z]);
    if lon < 0
        lon = lon+2*pi;
    end
    fprintf(fidOut, '#  %-8s  %13.4f %13.4f %13.4f  %8.4f %8.4f %6.1f\n', antenna(stat).name, antenna(stat).x, antenna(stat).y, antenna(stat).z, lat*180/pi, lon*180/pi, h);
end

% write comment line
fprintf(fidOut, '#\n#\n#\n#%s\n','                  Scan  YYYY.MM.DD-hh:mm:ss.s  Station   Azi(°)    Elev(°)  P(mbar) T[°C]  SlantPathDel(sec)');

% write data lines to file
for k = 1 : size(azel_data{1},1)
    
    % make proper format of delay
    totDelStr           = sprintf('%16.8d', ion(k));
%     totDelStr(end-2)    = []; % Why removing the "-" ???????
    
    fprintf(fidOut,  'O  %-16s%5.0f %4.0f.%02.0f.%02.0f-%02.0f:%02.0f:%04.1f  %-8s  %9.5f %8.5f  %6.1f %5.1f  %s\n', ['$',sessionName], azel_data{1}(k), years2write(k,1), months2write(k,1), days2write(k,1), azel_data{5}(k), azel_data{6}(k), azel_data{7}(k), azel_data{8}{k}, azel_data{9}(k)*180/pi, azel_data{10}(k)*180/pi, azel_data{13}(k), azel_data{12}(k), totDelStr);
end

% close file
fclose all;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% return error message of 0
errorMessage = 0;
varargout(1) = {errorMessage};
