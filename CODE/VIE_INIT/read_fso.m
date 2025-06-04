% ************************************************************************
%   Description:
%   Reads the fso file and saves the parameters in the struct orbit_data.
%
%
%   Input:										
%     satOrbitFilePathName      name of orbit file (FSO file)
%
% 
%   Output:
%     orbit_data                struct orbit data
% 
%   External calls: 
%	readstd.m
%       
%       
%   Coded for VieVS: 
%   17 Dec 2024 by Helene Wolf
%
%   Revision: 
%
%
% ************************************************************************
function [orbit_data] = read_fso(satOrbitFilePathName)

    [coeff,t0,hstep,nint,nq,nsat,tbound,tosc,oscele,satnum,narc,descr,source] = readstd(satOrbitFilePathName);
    
    orbit_data = struct('file_type', 'fso', 'file_pathname', satOrbitFilePathName, 'coeff', coeff, 't0', t0, 'hstep', hstep, 'nint', nint, 'nq', nq, 'nsat', nsat, 'tbound', tbound, ...
        'tosc', tosc, 'oscele', oscele, 'satnum', satnum, 'narc', narc, 'descr', descr, 'source', source);
end
