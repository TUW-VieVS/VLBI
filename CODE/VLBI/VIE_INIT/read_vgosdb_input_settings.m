% #########################################################################
% #     read_vgosdb_input_settings
% #########################################################################
%
% DESCRIPTION
% 	This file reads the vievs_input_settings.txt file from the /WORK/ folder
% 	and looks for "institute" and "frequencyband" path to file: ptf = 'vgosdb_input_settings.txt';
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
%	- ptf (path to file)
% OUTPUT:
%	- in (institute name)
%	- fb (frequency band name)
%   - wrapper_k  (wrapper tag, "all" or "ngs")
%   - wrapper_v  (wrapper version number)
%
% CHANGES:
%

function [ in, fb, wrapper_k, wrapper_v ] = read_vgosdb_input_settings( ptf )

in_tag          = 'institute:';         % => in
fb_tag          = 'frequency_band:';    % => in
wrapper_k_tag   = 'wrapper_k:';         % => wrapper_k
wrapper_v_tag   = 'wrapper_version:';   % => wrapper_v


in = '';
fb = '';
wrapper_k = '';
wrapper_v = '';

% ##### Open file #####
fid = fopen(ptf,'r');

if fid<0
    fprintf('File vievs_input_settings.txt does not exist\n')
else
    fprintf('Reading: %s\n', ptf)
    
    % ##### Loop over all lines #####
    while ~feof(fid)
	
        % #### Read one line ####
        str = fgetl(fid);
  
        if ~isempty(str)
            % parse line
            temp_str = textscan(str, '%s', 'CommentStyle', '%'); 
            if size(temp_str{1}, 1) == 2 % not a comment

                switch(temp_str{1}{1})

                    case in_tag
                        in = temp_str{1}{2};

                    case fb_tag
                        fb = temp_str{1}{2};

                    case wrapper_k_tag
                        wrapper_k = temp_str{1}{2};

                    case wrapper_v_tag
                        wrapper_v = temp_str{1}{2};

                end
            end
            
        end % if ~isempty(str)
        
    end % while ~feof(fid)
    
    
    
    
    % Close file
    fclose(fid);

    
end




