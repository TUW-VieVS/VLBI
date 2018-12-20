% ************************************************************************
%   Description:
%   Creates the parameter refname, refantbr and antactiv without reduced
%   antennas
%
%   Input:
%     anames_fo          antennas which will be session-wise reduced
%     refname_all        names of antennas
%     refantbr_all       array with all breaks and respective coordinates
%                        for all stations
%     antactiv_all       matrix with 1/0; columns are sessions; rows are
%                        antennas. 1 - antenna observed in this session, 0 - antenna didn't take
%                        part in this session
%
%   Output:                
%     refname            names of antennas in global adjustment
%     refantbr           array with breaks and respective coordinates
%                        for stations in global adjustment
%     antactiv           matrix with 1/0; columns are sessions; rows are
%                        antennas. 1 - antenna observed in this session, 0 - antenna didn't take
%                        part in this session. Only for antennas in global adjustment
%       
%   Coded for VieVS: 
%   21 Jul 2011 by Hana Spicakova
%
%   Revision: 
%%




function [refname, refantbr, antactiv]=refname_wored(anames_fo,refname_all,refantbr_all,antactiv_all)

    idcoor_fo=[];  % indices of stations to be reduced (wrt refname)
    nfo=size(anames_fo,1);
    for i=1:nfo
        [inant,xx]=find(strcmp(cellstr(anames_fo(i,1:8)),cellstr(refname_all)) == 1);
        idcoor_fo=[idcoor_fo; inant];
    end
    idcoor_fo=unique(idcoor_fo);
    
    % stations which will be estimated in the global adjustment (coordinates)
    stdf=setdiff([1:size(refname_all,1)],idcoor_fo);
    refname=refname_all(stdf,:);
    refantbr=refantbr_all(stdf);
    antactiv=antactiv_all(stdf,:);
    
    clear stdf
