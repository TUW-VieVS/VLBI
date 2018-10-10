% ************************************************************************
%   Description:
%   function to exclude  models of right ascension and declination as sourcewise 
%
%   Reference: 
%
%   Input:	
%       'opt'             structure array       (for info. /DOC/opt.doc)
%       'A_session'       structure array       design matrix of the real-observation equations
%       'num_inter_sou'  (1, num. sources)      number of estimation intervals for source coord.
%       'ns'             (1,1)                  number of sources
%
%   Output:
%       'A_session'     structure array     design matrix of the real-observation equations
%       'nso'           structure array     the number of source coordinate offsets for each source
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   30 Aug 2012 by Claudia Tierno Ros
%
%   Revision: 
%   19 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with
%   NNR condition
% ************************************************************************

function [A_session] = delsource_scan_A(opt,A_session,num_inter_sou,ns)
if opt.est_sourceNNR==0 
    for isou=ns:-1:1
        if   opt.source(isou).rade_inc == 0          

             if isou==1
                A_session(11).sm(:,1:num_inter_sou(1)+1) =[];
                A_session(12).sm(:,1:num_inter_sou(1)+1) =[];
             else
                A_session(11).sm(:,sum(num_inter_sou(1:isou-1))+isou:sum(num_inter_sou(1:isou))+isou) = [];
                A_session(12).sm(:,sum(num_inter_sou(1:isou-1))+isou:sum(num_inter_sou(1:isou))+isou) = [];
             end
             nso(isou).sources=0;

        end
    end
end
