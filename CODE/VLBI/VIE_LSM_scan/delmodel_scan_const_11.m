% ************************************************************************
%   Description:
%   function to exclude models from station coordinates pwl design matrix, weight matrix,
%   and o-c vector of constraints
%
%   Reference: 
%
%   Input:	
%       'HX'               structure array     design matrix of the pseudo-observation equations (constraints)
%       'PhX'              structure array     weight matrix of the pseudo-observation equations 
%       'ochX'             structure array     o-c vector of constraints (zero vector)
%       'opt'              structure array     (for info. /DOC/opt.doc)
%
%   Output:
%       'HX'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'PhX'         structure array     weight matrix of the pseudo-observation equations 
%       'ochX'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%
%   Revision: 
%
% ************************************************************************


function  [H11,Ph11,och11,H12,Ph12,och12 ]=delmodel_scan_const_11(ns,opt,num_inter_sou,sumso,H11,Ph11,och11,H12,Ph12,och12)

if opt.est_sourceNNR==0
    id1=[];
    id2=[];
    for iso=1 : ns 
        if opt.source(iso).rade_inc == 0    

            if iso==1
                id1=sumso.sources(iso)+1:sumso.sources(iso+1);
                id2=[1 sum(num_inter_sou(1:iso))];
            else
                id1 = [id1 sumso.sources(iso)+1:sumso.sources(iso+1)];
                id2 = [id2 sum(num_inter_sou(1:iso-1))+1:sum(num_inter_sou(1:iso))];
                
            end
        end
    end
    
    H11(:,id1) = [];
    H11(id2,:) = [];
    Ph11(id2,:) = [];
    Ph11(:,id2) = [];
    och11(id2)=[];

    H12(:,id1) = [];
    H12(id2,:) = [];
    Ph12(id2,:) = [];
    Ph12(:,id2) = [];
    och12(id2)=[];

end