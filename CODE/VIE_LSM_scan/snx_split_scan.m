% ************************************************************************
%   Description:
%   function which finds columns for reduction and columns which will be written
%   into the SINEX file
%
%   Reference: 
%
%   Input:										
%       x_           information about the columns
%       outsnx       1/0
%
%   Output:
%       col_est      vector with the parameters going into SINEX file
%       col_red      vector with the parameters which will be reduced from
%                    the N_sinex
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   11 oct 2012 Claudia Tierno from 01 May 2011 by Hana Spicakova
% ************************************************************************ 


function [col_est] = snx_split_scan(A_session,opt)

A1  = size(A_session(1).sm,2);
A2  = size(A_session(2).sm,2);
A3  = size(A_session(3).sm,2);
A4  = size(A_session(4).sm,2);
A5  = size(A_session(5).sm,2);
A6  = size(A_session(6).sm,2);
A7  = size(A_session(7).sm,2);
A8  = size(A_session(8).sm,2);
A9  = size(A_session(9).sm,2);
A10 = size(A_session(10).sm,2);
A11 = size(A_session(11).sm,2);
A12 = size(A_session(12).sm,2);
A13 = size(A_session(13).sm,2);
A14 = size(A_session(14).sm,2);
A15 = size(A_session(15).sm,2);
At  = A1+ A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9 + A10 + A11 + A12 ;


% clock parameters will be always reduced, in this version they cannot
% be written into SINEX file


% station coordinates cannot be reduced in SINEX file in this version
col_est = [ At+1: At+A13+A14+A15 ];
       

if opt.outsnx.zwd ~= 0
   
    col_est = [col_est [A1+A2+1:A1+A2+A3]];
end

if opt.outsnx.tgr ~= 0
    
    col_est = [col_est [A1+A2+A3+1:A1+A2+A3+A4+A5] ];
end

if opt.outsnx.sou ~= 0
   
    col_est = [col_est [A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+1:At]];
end

if opt.outsnx.eop ~= 0
   
    col_est = [col_est [A1+A2+A3+A4+A5+1:A1+A2+A3+A4+A5+A6+A7+A8+A9+A10]];
end


