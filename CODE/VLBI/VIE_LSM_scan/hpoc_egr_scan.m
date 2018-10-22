% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   egr
%
%   Reference: 
%
%   Input:	
%       'num_inter_egr'        (1,num. of antennas)     number of estimation intervals egr
%       'opt'                  structure array         (for info. /DOC/opt.doc)
%       'na'                   (1,1)                    number of antennas
%       'col_A5'               (1,1)                    number of columns of A_session(5)
%
%   Output:
%       'HX'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'PhX'         structure array     weight matrix of the pseudo-observation equations 
%       'ochX'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   25 July 2012 by Claudia Tierno Ros
%
%   Revision: 
%  
% ************************************************************************

function [H5,Ph5,och5] = hpoc_egr_scan(opt,na,col_A5,num_inter_egr)

Hrelegr = []; Phrelegr = [];
Habsegr = []; Phabsegr = [];
t_egr = sum(num_inter_egr); % total estimate interval of east gradients in a session   

for istat = 1:na 
    % FORMING THE RELATIVE CONSTRAINTS of east GRADIENTS
    coef_rel_egr = opt.stat(istat).coef_rel_egr;
    int_rel_egr = opt.stat(istat).int_egr;
   H_rel_egr(t_egr,num_inter_egr(istat)+1)= 0;
     P_rel_egr(t_egr,num_inter_egr(istat))= 0;
        sumrelegr = 0;

        sumrelegr=sum(num_inter_egr(1:istat-1));
        for inter =  1:num_inter_egr(istat) 
            H_rel_egr(sumrelegr + inter, inter) = +1;  % design matrix for the east gradients pseudo observtaions as constraints
            H_rel_egr(sumrelegr + inter, inter+1) = -1; % design matrix for the east gradients pseudo observtaions as constraints
            P_rel_egr(sumrelegr + inter, inter) = 1./coef_rel_egr^2; % [1/cm^2] weight matrix of the design matrix for the North gradients constraints
        end     
    Hrelegr = horzcat(Hrelegr,H_rel_egr); % Concatenating
    Phrelegr = horzcat(Phrelegr,P_rel_egr); % Concatenating

    % FORMING THE ABSOLUTE CONSTRAINTS of east GRADIENTS
    coef_abs_egr = opt.stat(istat).coef_abs_egr;
     
        for inter=1:(t_egr+na)
            Habsegr(inter, inter) = +1;  % design matrix for the east gradients pseudo observtaions as constraints
            Phabsegr (inter,inter)=1./coef_abs_egr^2; % [1/cm^2] weight matrix of the design matrix for the North gradients constraints
        end
         
    
    clear H_rel_egr P_rel_egr 
end

% FORMING THE O-C VECTOR FOR THE RELATIVE & ABSOLUTE CONSTRAINTS
och_rel(size(Hrelegr,1),1) = 0;  % O-C vector for the relative constraints of  east gradients
och_abs(size(Habsegr,1),1) = 0;  % O-C vector for the absolute constraints of  east gradients

%..........................................................................

% ELIMINATE PARAMETERS ACCORDING TO opt

for istat=na:-1:1

   if opt.stat(istat).egr_inc == 0
        
        if opt.constr_rel_egr == 1 && opt.constr_abs_egr == 0
            Hrelegr(:,sum(num_inter_egr(1:istat-1))+istat:sum(num_inter_egr(1:istat))+istat) = [];
            Hrelegr(sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat)),:) = [];
            
            Phrelegr(:,sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat))) = [];
            Phrelegr(sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat)),:) = [];
            
            och_rel(sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat))) = [];
        end

        if opt.constr_rel_egr == 0 && opt.constr_abs_egr == 1
            Habsegr(:,sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat))) = [];
            Habsegr(sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat)),:) = [];
            
            Phabsegr(:,sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat))) = [];
            Phabsegr(sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat)),:) = [];
            
            och_abs(sum(num_inter_egr(1:istat-1))+istat:sum(num_inter_egr(1:istat))+istat) = [];
        end
        
        if opt.constr_rel_egr == 1 && opt.constr_abs_egr == 1            
            Hrelegr(:,sum(num_inter_egr(1:istat-1))+istat:sum(num_inter_egr(1:istat))+istat) = [];
            Hrelegr(sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat)),:) = [];
           
            Phrelegr(:,sum(num_inter_egr(1:istat-1)):sum(num_inter_egr(1:istat))) = [];
            Phrelegr(sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat)),:) = [];
           
            och_rel(sum(num_inter_egr(1:istat-1)):sum(num_inter_egr(1:istat))) = [];

            Habsegr(:,sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat))) = [];
            Habsegr(sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat)),:) = [];
            
            Phabsegr(:,sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat))) = [];
            Phabsegr(sum(num_inter_egr(1:istat-1))+1:sum(num_inter_egr(1:istat)),:) = [];
            
            och_abs(sum(num_inter_egr(1:istat-1))+istat:sum(num_inter_egr(1:istat))+istat) = [];
        end
       
    end



end


%..........................................................................
% TAKE RELATIVE OR ABSOLUTE ACCORDING TO opt
if opt.constr_rel_egr==1 && opt.constr_abs_egr==0
   
    H5=Hrelegr;
    Ph5=Phrelegr;
    och5=och_rel;
   
end

if opt.constr_rel_egr==0 && opt.constr_abs_egr==1
  
    H5=Habsegr;
    Ph5=Phabsegr;
    och5=och_abs;
   
end

if opt.constr_rel_egr==1 && opt.constr_abs_egr==1
   
    H5=vertcat(Hrelegr,Habsegr);
    Ph5=blkdiag(Phrelegr,Phabsegr);
    och5=vertcat(och_rel,och_abs);
   
end

if opt.constr_rel_egr==0 && opt.constr_abs_egr==0
   
    if col_A5==0
       H5=[];
    else H5(1,col_A5)=0;
    end    
    Ph5=[];
    och5=[];
   
end

