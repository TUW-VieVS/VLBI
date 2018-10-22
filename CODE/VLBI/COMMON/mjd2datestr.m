% #########################################################################
% #     mjd2datestr
% #########################################################################
%
% DESCRIPTION
%   This function coverts a MJD vaue to the according date string.
%
%
% CREATED  
%   2015-08-28     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%
%
% INPUT
% - mjd              - Modified Julan Date
%
%
% OUTPUT
% - date_str         - Date string
%
% CHANGES:
%
function [date_str] = mjd2datestr(mjd, varargin)

    switch(nargin)
        
        case 1 % Default: Rounf to integer seconds
            flag_noRoundSec = 0;
            
        case 2 % Do not round seconds
            flag_noRoundSec = 1;
            
        otherwise
            fprintf('ERROR: Invalid number of inout arguments\n');
            return;
            
    end

    [y,m,d,h,min,sec] = mjd2date(mjd, flag_noRoundSec);

    if flag_noRoundSec
        fprintf('sec = %1.10f\n', sec);
    end

    date_str = datestr([y,m,d,h,min,sec]);
    % fprintf('sec = %10.8f\n', sec);
end

