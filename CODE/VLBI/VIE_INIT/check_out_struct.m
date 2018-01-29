% #########################################################################
% #     getinstName
% #########################################################################
%
% DESCRIPTION
% 	check_out_struct is a function to check and prepare the 
% 	out_struct for further processing. This involves:
% 	A.) Check if required files are available
% 	B.) Check version number of required files
% 	C.) Check institute provider of the required files
% 	D.) Distinguishes between different frequency bands
% 	The following files are used in VieVS and the reference to the VieVS
% 	out_struct is given:
% 	CrossReference/ObsCrossRef.nc --> out_struct.CrossReference.ObsCrossRef
% 	CrossReference/SourceCrossRef.nc --> out_struct.CrossReference.SourceCrossRef
% 	CrossReference/StationCrossRef.nc --> out_struct.CrossReference.StationCrossRef
% 	ObsEdit/GroupDelayFull_bX.nc --> out_struct.ObsEdit.GroupDelayFull_bX
% 	ObsEdit/Edit.nc --> out_struct.ObsEdit.Edit
% 	Observables/GroupDelay_bX.nc --> 
% 	ObsDerived/Cal_SlantPathIonoGroup_bX.nc --> out_struct.ObsDerived.Cal_SlantPathIonoGroup_bX
% 	Scan/TimeUTC.nc
% 	<Station>/TimeUTC.nc
% 	<Station>/Met.nc
% 	<Station>/Cal_Cable.nc
%   
%
% CREATED  
%   2017-11-14     Jakob Gruber
%
% REFERENCES
%
%
% COUPLING
%
%
% INPUT:
%	- out_struct (matlab struct of the vgos data base)
%       - in (string institute name)
% OUTPUT:
%	- out_struct_checked (matlab struct of the checked vgos data base)
%
% CHANGES:
%
function [ out_struct_checked ] = check_out_struct( out_struct, in )

verbosity = 0;

%% INPUT
dn{1} = 'ObsEdit';files{1} = {'GroupDelayFull','Edit'};
dn{2} = 'CrossReference';files{2} = {'ObsCrossRef','SourceCrossRef','StationCrossRef'};
dn{3} = 'Observables'; files{3} = {'GroupDelay'};
dn{4} = 'ObsDerived'; files{4} = {'Cal_SlantPathIonoGroup'};
dn{5} = 'Scan'; files{5} = {'TimeUTC'};
dn_stat = 'stat'; files_stat = {'TimeUTC','Met','Cal_Cable'};


%% static

fnames = fieldnames(out_struct);

id_V = '_V';        % identifier for version number
id_inst = '_i';     % identifier for institute name
id_fband = '_b';   % identifier for frequency band

%% loop through whole structure
for ifn = 1:numel(fnames) 
    %% loop trough specified folders dn
    for idn = 1:length(dn)
        %% lets find specified folder in the whole data base
        if strcmp(fnames{ifn},dn{idn})
            fn = fieldnames(out_struct.(fnames{ifn})); % get fieldnames of sub structure of db
            %% loop through desired files of current folder
            for i1 = 1:length(files{idn})
                
                %% find substructures which correspond to the desired files
                ifi = find(strncmp(fn,files{idn}{i1},length(files{idn}{i1})));
                tmp = fn(ifi); % found filenames of desired substructures
                %% Prealllocation of version number vector, institute name and frequency band name cell
                num_curr_ver = -99*ones(length(tmp),1);
                name_curr_inst = cell(length(tmp),1);
                name_curr_freq = cell(length(tmp),1);
                %% Get version number, institute name and frequency band name of each file name
                % if version number does not exist --> version number is set to -1
                % if insitute name does not exist --> corresponding element in the cell array is empty
                % if frequency band name does not exist --> corresponding element in the cell array is empty
                for i = 1:length(tmp)
                    num_curr_ver(i)  = getVerNr( tmp{i},id_V );
                    name_curr_inst{i} = getInstName( tmp{i},id_inst );
                    name_curr_freq{i} = getFreqBand( tmp{i},id_fband );
                end
                
                if isempty([name_curr_freq{:}])
                    %% not needed to distinguish between frequency bands                    
                    IFB = {'default'}; % set default value for frequency band cell if no frequency band exist
                    nloop = 1;
                else
                    %% Loop over several frequency bands
                    IFB = unique(name_curr_freq);
                    nloop = length(IFB);
                end
                for ifband = 1:nloop
                    fbcurr = IFB{ifband};
                    
                    %% Find frequency band name
                    ifbcurr    = find(strcmp(fbcurr,name_curr_freq));
                    %% Find instite name
                    iinst = find(strcmp(in,name_curr_inst));
                    %% Find version name
                    ivern = find(num_curr_ver~=-1);                    
                    %% Special case
                    if ~isempty(ifbcurr) && ~isempty(ivern)
                        if isempty(intersect(ifbcurr,ivern))
                            ivern = [];
                        end                        
                    end
                    %% Check all possible cases: fband = 0/1 , inst = 0/1 , vernr = 0/1
                    if isempty(ifbcurr)    % fband = 0
                        if isempty(iinst)   % fband = 0, inst = 0
                            if isempty(ivern)
                                %% fband = 0, inst = 0, version = 0
                                tmp_out = fn{ifi};
                            else
                                %% fband = 0, inst = 0, version = 1
                                [~,imaxver] = max(num_curr_ver); % find maximum version number
                                tmp_out = fn{ifi(imaxver)}; % store file                                
                            end
                        else % fband = 0, inst = 1
                            if isempty(ivern)
                                %% fband = 0, inst = 1, version = 0
                                iii = fn{ifi(iinst)};
                                tmp_out = iii; % store file                                
                            else
                                %% fband = 0, inst = 1, version = 1
                                ii = intersect(ivern,iinst);
                                [~,imaxver] = max(num_curr_ver(ii)); % find maximum version number
                                tmp_out = fn{ifi(ii(imaxver))}; % store file
                            end
                        end
                    else    % fband = 1
                        if isempty(iinst)  % fband = 1, inst = 0
                            if isempty(ivern)
                                %% fband = 1, inst = 0, version = 0
                                tmp_out = fn{ifi(ifbcurr)}; % store file
                            else
                                %% % fband = 1, inst = 0, version = 1
                                ii = intersect(ifbcurr,ivern);
                                [~,imaxver] = max(num_curr_ver(ii)); % find maximum version number
                                tmp_out = fn{ifi(ii(imaxver))}; % store file
                            end
                        else % fband = 1, inst = 1
                            if isempty(ivern)
                                %% fband = 1, inst = 1, version = 0
                                ii = intersect(ifbcurr,iinst);
                                tmp_out = fn{ifi(ii)}; % store file
                            else
                                %% fband = 1, inst = 1, version = 1
                                ii = intersect(ifbcurr,iinst);
                                iii = intersect(ii,ivern);
                                [~,imaxver] = max(num_curr_ver(iii)); % find maximum version number
                                tmp_out = fn{ifi(iii(imaxver))}; % store file
                            end
                        end
                    end
                    if 1 == 1
                        fprintf('Load file %s \n',tmp_out)
                    end
                    if isempty(ifbcurr)
                        out_struct.(dn{idn}).(files{idn}{i1}) = out_struct.(dn{idn}).(tmp_out); % put temporarly stored file to data
                    else
                        if strcmp(fbcurr,'X')
                            out_struct.(dn{idn}).([files{idn}{i1},'_bX']) = out_struct.(dn{idn}).(tmp_out); % put temporarly stored file to data
                        end
                        if strcmp(fbcurr,'S')
                            out_struct.(dn{idn}).([files{idn}{i1},'_bS']) = out_struct.(dn{idn}).(tmp_out); % put temporarly stored file to data
                        end
                    end
                end
            end                
        end
        %% Check Station Folders
        if strcmp(fnames{ifn},dn_stat)
            %% Loop through all stations
            for ist = 1:numel(out_struct.(fnames{ifn}))
                fn = fieldnames(out_struct.(fnames{ifn})(ist));
                %% loop through desired files of current folder
                for i1 = 1:length(files_stat)
                    %% find substructures which correspond to the desired files
                    ifi = find(strncmp(fn,files_stat{i1},length(files_stat{i1})));
                    tmp = fn(ifi); % found filenames of desired substructures
                    %% check if file exists
                    if isempty(tmp)
                        fprintf('File missing: %s\n',files_stat{i1})
                        continue;
                    end
                    %% Prealllocation of version number vector, institute name and frequency band name cell
                    num_curr_ver = -99*ones(length(tmp),1);
                    name_curr_inst = cell(length(tmp),1);
                    name_curr_freq = cell(length(tmp),1);
                    %% Get version number, institute name and frequency band name of each file name
                    % if version number does not exist --> version number is set to -1
                    % if insitute name does not exist --> corresponding element in the cell array is empty
                    % if frequency band name does not exist --> corresponding element in the cell array is empty
                    for i = 1:length(tmp)
                        num_curr_ver(i)  = getVerNr( tmp{i},id_V );
                        name_curr_inst{i} = getInstName( tmp{i},id_inst );
                        name_curr_freq{i} = getFreqBand( tmp{i},id_fband );
                    end
                    
                    if isempty([name_curr_freq{:}])
                        %% not needed to distinguish between frequency bands
                        IFB = {'default'}; % set default value for frequency band cell if no frequency band exist
                        nloop = 1;
                    else
                        %% Loop over several frequency bands
                        IFB = unique(name_curr_freq);
                        nloop = length(IFB);
                    end
                    for ifband = 1:nloop
                        fbcurr = IFB{ifband};                        
                        %% Find frequency band name
                        ifbcurr    = find(strcmp(fbcurr,name_curr_freq));
                        %% Find instite name
                        iinst = find(strcmp(in,name_curr_inst));
                        %% Find version name
                        ivern = find(num_curr_ver~=-1);
                        %% Special case
                        if ~isempty(ifbcurr) && ~isempty(ivern)
                            if isempty(intersect(ifbcurr,ivern))
                                ivern = [];
                            end                        
                        end
                        %% Check all possible cases: fband = 0/1 , inst = 0/1 , vernr = 0/1
                        if isempty(ifbcurr)    % fband = 0
                            if isempty(iinst)   % fband = 0, inst = 0
                                if isempty(ivern)
                                    %% fband = 0, inst = 0, version = 0
                                    tmp_out = fn{ifi};
                                else
                                    %% fband = 0, inst = 0, version = 1
                                    [~,imaxver] = max(num_curr_ver); % find maximum version number
                                    tmp_out = fn{ifi(imaxver)}; % store file                                
                                end
                            else % fband = 0, inst = 1
                                if isempty(ivern)
                                    %% fband = 0, inst = 1, version = 0
                                    iii = fn{ifi(iinst)};
                                    tmp_out = iii; % store file                                
                                else
                                    %% fband = 0, inst = 1, version = 1
                                    ii = intersect(ivern,iinst);
                                    [~,imaxver] = max(num_curr_ver(ii)); % find maximum version number
                                    tmp_out = fn{ifi(ii(imaxver))}; % store file
                                end
                            end
                        else    % fband = 1
                            if isempty(iinst)  % fband = 1, inst = 0
                                if isempty(ivern)
                                    %% fband = 1, inst = 0, version = 0
                                    tmp_out = fn{ifi(ifbcurr)}; % store file
                                else
                                    %% % fband = 1, inst = 0, version = 1
                                    ii = intersect(ifbcurr,ivern);
                                    [~,imaxver] = max(num_curr_ver(ii)); % find maximum version number
                                    tmp_out = fn{ifi(ii(imaxver))}; % store file
                                end
                            else % fband = 1, inst = 1
                                if isempty(ivern)
                                    %% fband = 1, inst = 1, version = 0
                                    ii = intersect(ifbcurr,iinst);
                                    tmp_out = fn{ifi(ii)}; % store file
                                else
                                    %% fband = 1, inst = 1, version = 1
                                    ii = intersect(ifbcurr,iinst);
                                    iii = intersect(ii,ivern);
                                    [~,imaxver] = max(num_curr_ver(iii)); % find maximum version number
                                    tmp_out = fn{ifi(iii(imaxver))}; % store file
                                end
                            end
                        end
                        if verbosity == 1
                            fprintf('Load file %s \n',tmp_out)
                        end
                        %% In case data is empty do not update vievs internal struct
                        if ~isempty(out_struct.(dn_stat)(ist).(tmp_out)) 
                            out_struct.(dn_stat)(ist).(files_stat{i1}) = out_struct.(dn_stat)(ist).(tmp_out); % put temporarly stored file to data
                        end 
                    end
                    
                end
            end
        end
    end
end

out_struct_checked = out_struct;


end

