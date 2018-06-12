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
%
% CHANGES:
%

function [ in,fb ] = read_vgosdb_input_settings( ptf )

inid = 'institute:';
fbid = 'frequencyband:';

in = [];
fb = [];

fid = fopen(ptf);
if fid<0
    fprintf('File vievs_input_settings.txt does not exist\n')
else
    fprintf('File vievs_input_settings.txt found\n')
    tline = fgetl(fid);
    while ischar(tline)
        ic = strfind(tline,'%'); % find comments
        if ~isempty(ic)
            tline = tline(1:ic-1);
        end
        is = strfind(tline,' ');
        if ~isempty(is)
            %         for iis = 1:length(is)
            %             ec{iis} = '';
            %         end
            tline(is) = '';
        end
        iin = strfind(tline,inid); % find institute
        ifb = strfind(tline,fbid); % find frequency band
        if ~isempty(iin)
            in = tline(iin+length(inid):end);
        end
        if ~isempty(ifb)
            fb = tline(ifb+length(fbid):end);
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    
    if isempty(in)
        fprintf('No institute is specified\n')
    else
        fprintf('Institute specified: %s\n',in)
    end
    if isempty(fb)
        fprintf('No frequencyband is specified\n')
    else
        fprintf('Frequencyband specified: %s\n',fb)
    end    
end


end

