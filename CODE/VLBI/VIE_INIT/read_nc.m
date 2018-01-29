%% This function reads netCDF files
%Status quo:
%       Reads out nc files from a VLBI session and passes a struct
%       containing structs of each station as well as one struct for 'obs'
%       and one for 'head'.
%       That is due to the fact that each session must have a head.nc as
%       well as an 'OBS' folder in the specified directory. 
%
%
%Created:
%       Matthias Madzak, March 11th 2010
%
% in:
%       head.nc     (string)    full path of head.nc file (or only path to
%                               head/Head.nc file)
%
% out:
%       out_struct  (1x1 struct)containing station and observation folders 
%                               and the head struct
%                               --eg---------------------------------------
%                               out_struct.stat(1) eg wettzell     (1x1 struct)
%                               out_struct.stat(2) eg westford     (1x1 struct)
%                               out_struct.stat(...)               (1x1 struct)
%
%       nc_info     (1xnstat+1) struct containing information about the
%                               session
%                               --eg---------------------------------------
%                               nc_info(1).foldername       ('kokee')
%                               nc_info(1).file_list        ('AzEl'; 'CableCal'; 'Cal_kAxisOffSet'; 'Met'; 'Part_kAxisOffSet')
%                               nc_info(1).variable_list    (nfilesx20 cell) ->
%                                       for each file in folder one row, in cols: variables of file(row) (if fewer than 20 variables: rest of cols are empty)
%                                       to get number of variables from 
%                                       folder x and file y:
%                                       length(find(~cellfun('isempty', {nc_info(1,x).variable_list{y,:}})))
%
%
%   2016-06-16, M. Madzak: Errors when loading intensive sesssions: - Error, when variable name ("varname_new") contained "+" => Replaced with "_" in "varname_new"
%   2017-11-13, J. Gruber: Now the version number in the header file (Head.nc) which is stored in the top level of the vgosDb of a session is considered. By default the highest version will be used for processing.
%   2017-11-14, J. Gruber: Stations with spaces in their names can be processed

function [out_struct, nc_info]=read_nc(headNcFile)
% define head-filename, !! should NOT be changed !!
%headname='head.nc';

% define folders which are read IF THEY EXIST (in addition to station
% folders)
folders2read={'Apriori', 'CrossReference', 'History', 'ObsCalTheo', 'ObsDerived', ...
    'ObsEdit', 'Observables', 'ObsPart', 'ObsTheoretical', 'Scan',...
    'Session', 'Solve'};

%make sure name ends with .nc
% if length(name)>=3
%     if ~strcmp(name(length(name)-2:length(name)), '.nc')
%         name=[name, '.nc'];
%     end
% else
%     name=[name, '.nc'];
% end

%% CHECKS / PREs
% get directory from full path -> directory
directory=strrep(strrep(headNcFile,'Head.nc',''), 'head.nc', '');
% Take Head.nc file with highest Version number
top_level = dir(directory);
for i_toplev = 1:length(top_level)
    tmp(i_toplev) = strncmp(top_level(i_toplev).name,'Head',4);
end
i_head_files = find(tmp);
Head_files={top_level(i_head_files).name};
id_ver = '_V';
for i_toplev_head = 1:length(Head_files)
    curr_headFile = Head_files{i_toplev_head};
    Head_version(i_toplev_head) = getVerNr( curr_headFile,id_ver );  
end
[~,i_max_ver]=max(Head_version);
headNcFile=dir([directory,top_level(i_head_files(i_max_ver)).name]);
headNcFile=[directory,'/',headNcFile(1).name];

% get all files from directory
dirContentOfSession=dir(directory);

% define mode to read nc-file
mode='NC_NOWRITE';

if isempty(dirContentOfSession)
    fprintf('No files found in %s\nReturning from ''read.nc''', directory)
    return
end

% check if there is an OBS folder
if ~isdir([directory, 'Observables'])
    fprintf('Folder ''Observables'' does not exist in ''%s''\nReturning from read_nc.m\n', directory)
    keyboard;
    return
end

%% read head.nc
if ~exist(headNcFile, 'file')
    fprintf('%s was not found!\nReturning from read_nc.m\n', headNcFile)
    return
end

% create struct
head=struct;

% open nc-file
ncid=netcdf.open(headNcFile, mode);

% get # of dimensions, variables, global attributes and index of
% unlimited dimension
[ndims,nvars,ngatts,unlimdimid]=netcdf.inq(ncid);

% PRINT - commented!
%fprintf('\nIn this file there are %1d dimensions, %1d variables and %1d global Attributes. \nThe dimension with unlimited length has ID: %1.d\n\n', ndims, nvars, ngatts, unlimdimid)

% for all variables in this .nc file make one field in struct array and
% load data into it
for varid=0:nvars-1
         
    % get infos from .nc file for 'running' varid
    [varname, xtype, varDimIDs, varAtts]=netcdf.inqVar(ncid, varid);
        
    % replace bad characters by '_'
    varname_new=strrep(varname, ':', '_');
    varname_new=strrep(varname_new, '@', '_');
        
    % read out data from variable
    eval(['head.', varname_new, '.val=netcdf.getVar(ncid, varid);']);
    
    % PRINT - commented!
    %fprintf('          Variable ''%11s'' ----> %10s.%s\n', varname, 'head', varname_new);
            
end

%close .nc file
netcdf.close(ncid);  


%% Create Output (1/2)

 

% get num of stations=num of folders
nStat=size(head.StationList.val',1);

% Check stationlist for spaces in strings
for istat = 1:nStat
    stm = strtrim(head.StationList.val(:,istat)');
    is = strfind(stm,' ');
    if ~isempty(is)
        head.StationList.val(is,istat) = '_'; % replace space with underscore ' ' --> '_'
    end        
end

% nc_info is a struct of length numberOfStations+1 (OBS) -> so for every
% folder in the session directory one entry.
% STATION/'OBS' -> FILE -> VARIABLE LIST
% create nc_info struct (1of2 overgiven variables)
nc_info(nStat+1)=struct('foldername', [], 'file_list', [], 'variable_list', []); 

%e.g. nc_info(1)
% .foldername: wettzell
% .file_list: 'Atm'; 'AzEl', ...
% .variable_list: 'Stub', 'createdBy', 'Session', 'Cable_1', ... -> for Atm
%                 'Stub', 'createdBy', 'Az',      'El',      ... -> for AzEl
%                                                                -> ... for
%                                                                each file
%                                                                one row of
%                                                                contained
%                                                                variables


%% check if all stations mentioned in head.nc have one folder in 'directory'
% and contain (at least) 'Met.nc' (just as an example)

% chek if all folders exist
for k=1:nStat  %for all stations
    if ~isdir([directory, strrep(head.StationList.val(1:8,k)',' ','')])       %all these folders should be contained in 'directory'
        fprintf('\nfolder ''%s'' not found in ''%s''\nReturning from read_nc.m\n', head.StationList.val(1:8,k)', directory)
        return
        
    % the following check (if Met.nc is in folder) can't be done because
    % '...\d08JAN02XA\KOKEE   \Met.nc' won't find anything
    
    %elseif ~exist([directory,  head.StatNames((k-1)*8+1:k*8), '\Met.nc'], 'file')
    %    fprintf('File ''Met.nc'' not found in ''%s''\nReturning from read_nc.m\n', [directory, head.StatNames((k-1)*8+1:k*8)])
    %    return
    end
end


%% reading .nc files from station folders

stat=struct;

for iStat=1:nStat %for all stations
    
    curStatName=strrep(head.StationList.val(1:8,iStat)',' ','');

    % create struct named like the station
    %eval([lower(curStatName), '=struct;']); % This is not good programming
    %practice! -> make variable name independent of stationname
    
    %eval([curStatName, '=struct;']);
        
    % get files in stationName-folder
    files=dir([directory, curStatName, '/*.nc']);
    
    % write stationname to nc_info (output) and create file_list cell and
    % variable_list cell. Even though nvars - which is num of cols in 
    % variable_list - is yet unknown. But max num of variables is (by now) 
    % 20 -> so create 20 columns.
    % HINT: get number of variables in one row of nc_info.variable_list:
    % length(find(~cellfun('isempty', {nc_info(1,2).variable_list{3,:}})))
    nc_info(iStat).foldername=curStatName;
    %nc_info(stat).foldername=curStatName;
    nc_info(iStat).file_list=cell(length(files), 1);
    nc_info(iStat).variable_list=cell(length(files), 20);
    
    for k=1:length(files)   %for all files (in folder of station name)
        % define file to be read
        curFile=files(k).name;
        
        struct_name=curFile(1 : strfind(curFile, '.')-1);
        struct_name=strrep(struct_name,'-','_');
        
%         eval([lower(curStatName), '.', struct_name, '=struct;']);
        if k==1
            stat(iStat).(struct_name)=struct;
        end
        
        % fill in struct/variable name
        nc_info(iStat).file_list{k, 1}=struct_name;
        
        % open nc-file
        ncid=netcdf.open([directory, curStatName, '/', curFile], mode);
        
        % get # of dimensions, variables, global attributes and index of
        % unlimited dimension
        [ndims,nvars,ngatts,unlimdimid]=netcdf.inq(ncid);
        
        % PRINT - commented!
        %fprintf('\nIn this file there are %1d dimensions, %1d variables and %1d global Attributes. \nThe dimension with unlimited length has ID: %1.d\n\n', ndims, nvars, ngatts, unlimdimid)
        
        % for all variables in this .nc file make one field in struct array and
        % load data into it
        for varid=0:nvars-1
            % get infos from .nc file for 'running' varid
            [varname, xtype, varDimIDs, varAtts]=netcdf.inqVar(ncid, varid);
            
            % replace bad characters by '_'
            varname_new=strrep(varname, ':', '_');
            varname_new=strrep(varname_new, '@', '_');
            varname_new=strrep(varname_new, '-', '_');
            
            % fill in variable name in struct
            nc_info(iStat).variable_list{k, varid+1}=varname_new;
        
            % read out data from variable
            stat(iStat).(struct_name).(varname_new).val=netcdf.getVar(ncid, varid);
   
%             eval([lower(curStatName), '.', struct_name, '.', varname_new, '=netcdf.getVar(ncid, varid);']);
            
            % PRINT - commented!
            %fprintf('          Variable ''%11s'' ----> %8s.%s.%s\n',
            %varname, lower(curStatName), struct_name, varname_new);
        
            % save information to variable_list (for every variable one entry)
            %variable_list(m).ncFile=curFile;          %.nc file the variable is coming from
            %variable_list(m).varnameInNcFile=varname; %variable name in .nc file (incl bad characters, like ':')
            %variable_list(m).name=[lower(curStatName), '.', struct_name, '.', varname_new];    %stored variable name (struct.field)
%             switch xtype %may be extended                   %save data type of variable
%                 case 2
%                     variable_list(m).type='char';
%                 case 3
%                     variable_list(m).type='int16';
%                 case 4
%                     variable_list(m).type='int32';
%                 case 5
%                     variable_list(m).type='single';
%                 case 6
%                     variable_list(m).type='double';
%                 otherwise
%                     variable_list(m).type='unknown';
%             end
            %m=m+1;
        end
        
        %close .nc file
        netcdf.close(ncid);  
    end
    
    % PRINT - commented!
    %fprintf('\n=========================================================================\n')
    %fprintf('=========================================================================\n\n')
    
end


%% reading .nc files from other but station folders
out_struct=struct;
for iFolder=1:length(folders2read)

    curFolderName=folders2read{iFolder};
    
    properFolderLogical=strcmpi({dirContentOfSession.name},curFolderName);
    
    % if the current folder does exist
    if sum(properFolderLogical)==1
        
        curFolderName=strrep(curFolderName,':','_');
        curFolderName=strrep(curFolderName,'@','_');
        curFolderName=strrep(curFolderName,'-','_');
        curFolderName=strrep(curFolderName,'+','_');
        
        % create struct named like the station
        out_struct.(curFolderName)=struct;
        % eval([lower(curStatName), '=struct;']);

        % get files in stationName-folder
        files=dir([directory, ...
            dirContentOfSession(properFolderLogical).name, '/*.nc']);

        % write stationname to nc_info (output) and create file_list cell and
        % variable_list cell. 
        nc_info(nStat+iFolder).foldername=curFolderName;
        %nc_info(nStat+1).foldername=curStatName;
        nc_info(nStat+iFolder).file_list=cell(length(files), 1);
        nc_info(nStat+iFolder).variable_list=cell(length(files), 20);

        for k=1:length(files)   %for all files (in folder of station name)
            % define file to be read
            curFile=files(k).name;

            % PRINT - commented!
            %fprintf('\nCurrent File: ''%s''\n', curFile)

            % create struct (with same name as .nc-file)
            %backsl_ind=findstr(curFile, '\');
            struct_name=curFile(1 : strfind(curFile, '.')-1);
            struct_name=strrep(struct_name,':','_');
            struct_name=strrep(struct_name,'@','_');
            struct_name=strrep(struct_name,'-','_');

            out_struct.(curFolderName).(struct_name)=struct;
        %     eval([lower(curStatName), '.', struct_name, '=struct;']);

            % fill in struct/variable name
            nc_info(nStat+iFolder).file_list{k, 1}=struct_name;

            % open nc-file
            ncid=netcdf.open([directory, ...
                dirContentOfSession(properFolderLogical).name, '/', curFile], mode);

            % get # of dimensions, variables, global attributes and index of
            % unlimited dimension
            [ndims,nvars,ngatts,unlimdimid]=netcdf.inq(ncid);

            % PRINT - commented!
            %fprintf('\nIn this file there are %1d dimensions, %1d variables and %1d global Attributes. \nThe dimension with unlimited length has ID: %1.d\n\n', ndims, nvars, ngatts, unlimdimid)

            % for all variables in this .nc file make one field in struct array and
            % load data into it
            for varid=0:nvars-1
                % get infos from .nc file for 'running' varid
                [varname, xtype, varDimIDs, varAtts]=netcdf.inqVar(ncid, varid);
                
                
                % replace bad characters by '_'
                varname_new=strrep(varname, ':', '_');
                varname_new=strrep(varname_new, '@', '_');
                varname_new=strrep(varname_new, '-', '_');
                varname_new=strrep(varname_new, '+', '_');

                % fill in variable name in struct
                nc_info(nStat+iFolder).variable_list{k, varid+1}=varname_new;

                % read out data from variable
                try
                    out_struct.(curFolderName).(struct_name).(varname_new).val=netcdf.getVar(ncid, varid);
                catch
                    happy_birthday
                   disp('test'); 
                end
  %         eval([lower(curStatName), '.', struct_name, '.', varname_new, '=netcdf.getVar(ncid, varid);']);
                
                % for all attributes of current variable
                for attnum=0:varAtts-1
                    % Get attribute name, given variable id.
                    attname = netcdf.inqAttName(ncid,varid,attnum);
                    attval  = netcdf.getAtt(ncid,varid,attname);
                    out_struct.(curFolderName).(struct_name).(varname_new).attr(attnum+1).name=...
                        attname;
                    out_struct.(curFolderName).(struct_name).(varname_new).attr(attnum+1).val=...
                        attval;
                end
                % PRINT - commented!
                %fprintf('          Variable ''%11s'' ----> %8s.%s.%s\n', varname, lower(curStatName), struct_name, varname_new);
            end

            %close .nc file
            netcdf.close(ncid);
        end   
    end
end

% PRINT - commented!
%fprintf('\n=========================================================================\n')
%fprintf('=========================================================================\n\n')


%% Saving all created structs (beside 'variable_list' which is itself 
% overgiven) into one struct that is overgiven.

% create output struct (= 2of2 overgiven variables)
% nfolders=length(nc_info); % this must be equal nStat+1 (1... observables)


for iStat=1:nStat
    out_struct.stat(iStat)=stat(iStat);
end
%out_struct.observables=observables;
%     eval(['out_struct.', nc_info(k).foldername, '=', nc_info(k).foldername, ';']);

out_struct(1).head=head;


