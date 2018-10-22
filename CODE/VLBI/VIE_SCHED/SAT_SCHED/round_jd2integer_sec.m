% #########################################################################
% #     round_jd2integer_sec
% #########################################################################
%
% DESCRIPTION
%   This function rounds JD epochs to integer seconds.
%   There is the option to round up or down.
%
% CREATED  
%   2016-08-02     Andreas Hellerschmied
%
% REFERENCES
%
% COUPLING
%
% INPUT
% - jd_in       - input D
% - method      - "round", "up" or "down" (optional)
%
%
% OUTPUT
% - jd_int_sec  - JD, rounded to integer seconds
%
% CHANGES:

function [jd_int_sec] = round_jd2integer_sec(jd_in, varargin)

% Init.:
jd2sec = 24*60*60;

switch(nargin)
    case 2
        method = varargin{1};
    case 1
        method = 'round';
    otherwise
        error('Invalid number of input arguments!');
end

% [y,m,d,h,min,sec] = mjd2date(jd_in-2.400000500000000e+006)

% Convert JD to sec:
sec_in = jd_in * jd2sec;

switch(method)
    case 'round'
        sec_int = round(sec_in);
    case 'down'
        sec_int = floor(sec_in);
    case 'up'
        sec_int = ceil(sec_in);
end

jd_int_sec = sec_int / jd2sec;

% [y,m,d,h,min,sec] = mjd2date(jd_int_sec-2.400000500000000e+006)

