% #########################################################################
% #     getTLEFileNames
% #########################################################################
%
% DESCRITPION
% Function to get the names of all files with the ending *.tle from the 
% directory defined in "tle_dir".
%
% AUTHOR 
%   Andreas Hellerschmied (2015.01.21)
%
% INPUT
%   path_tle_dir   Path to TLE file directory
%
% OUTPUT
%   tle_str         String containing all filenames
%
% CHANGES


function [tle_str] = getTLEFileNames(path_tle_dir)

    % --- get avilable TLE filenames ---
    tle_files_struct = dir([path_tle_dir ,'*.tle']);

    number_of_tle_files = length(tle_files_struct);
    tle_str = '';

    for i = 1 : (number_of_tle_files - 1)
        tle_str = [tle_str, tle_files_struct(i).name, '|'];
    end
    tle_str = [tle_str, tle_files_struct(number_of_tle_files).name];
    
end


