% #########################################################################
% #     jd2datestr
% #########################################################################
%
% DESCRIPTION
%   This function coverts a JD vaue to the according date string.
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
% - jd                  - Julan date
% - put_format          - Output format specification for the function "datestr" (string) - optional
%
%
% OUTPUT
% - date_str           - Date string
%
% CHANGES:
%
function [date_str] = jd2datestr(jd, varargin)

    switch(nargin)
        case 1
            out_format = 'yyyy-mm-dd HH:MM:SS';
        case 2
            out_format = varargin{1};
    end

    [y,m,d,h,min,sec] = mjd2date(jd-2.400000500000000e+006);
    date_str = datestr([y,m,d,h,min,sec], out_format);
end

