% ************************************************************************
%   Description:
%   function to rearrange the data from scan based to source based
%
%   Reference: 
%
%   Input:	
%       'opt'              structure array     (for info. /DOC/opt.doc)
%       'scan'             structure array     (for info. /DOC/scan.doc)
%       'iscan'            (1,1)                Number of the scan that is being analysed
%       'sources'          structure array     (for info. /DOC/sources.doc)
%       't'                (1,num. of scans)    minutes from midnight to the time of the scans
%       'mjd_scan'         (1,num. of scans)    mjd of the scans
%
%   Output:
%       'per_source'       structure array     (for info. /DOC/per_source.doc)
%       'T_source'         (1,2)                estimation epochs for clocks, zwd, ngr, egr, xyz, sou
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   July 2012 by CLaudia Tierno Ros
%
%   Revision: 
%  
% ************************************************************************

function [per_source,T_source] = sourcewisepar_scan(opt,scan,iscan,sources,t,mjd_scan)
 
    per_source.mjd=mjd_scan(iscan);
    per_source.minute=t(iscan);
    per_source.total=scan(iscan).nobs;
    per_source.iso=scan(iscan).iso;
    per_source.name=sources(per_source.iso).name;

    for j=1:length (scan(iscan).obs)
        per_source.ra(1,j)=scan(iscan).obs(j).psou(1);  % [cm/mas]
        per_source.de(1,j)=scan(iscan).obs(j).psou(2);  % [cm/mas]
    end
    
    per_source.no_int_rade=1;

    t1 = floor(t(iscan)/opt.sour_int_rade)*opt.sour_int_rade;
    t2 =ceil(t(iscan)/opt.sour_int_rade)*opt.sour_int_rade;
    T_source = [t1,t2]; % estimation intervals for source


