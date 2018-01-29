% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   ngr
%
%   Reference: 
%
%   Input:	
%       'num_inter_ngr'        (1,num. of antennas)     number of estimation intervals ngr
%       'opt'                  structure array         (for info. /DOC/opt.doc)
%       'na'                   (1,1)                    number of antennas
%       'col_A4'               (1,1)                    number of columns of A_session(4)
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


function [H4,Ph4,och4] = hpoc_ngr_scan(opt,na,col_A4,num_inter_ngr)

Hrelngr = []; Phrelngr = [];
Habsngr = []; Phabsngr = [];
t_ngr = sum(num_inter_ngr); % total estimate interval of north gradients in a session   

for istat = 1:na 
    % FORMING THE RELATIVE CONSTRAINTS of NORTH GRADIENTS
    coef_rel_ngr = opt.stat(istat).coef_rel_ngr;
    int_rel_ngr = opt.stat(istat).int_ngr;
    H_rel_ngr(t_ngr,num_inter_ngr(istat)+1)= 0;           
    P_rel_ngr(t_ngr,num_inter_ngr(istat))= 0;     
    

    sumrelngr=sum(num_inter_ngr(1:istat-1));
        for inter =1:num_inter_ngr(istat)                    
            H_rel_ngr(sumrelngr + inter, inter) = +1;  % design matrix for the North gradients pseudo observtaions as constraints
            H_rel_ngr(sumrelngr + inter, inter+1) = -1; % design matrix for the North gradients pseudo observtaions as constraints
            P_rel_ngr(sumrelngr + inter, inter) =  1./coef_rel_ngr^2; % [1/cm^2] weight matrix of the design matrix for the North gradients constraints
        end     
    Hrelngr = horzcat(Hrelngr,H_rel_ngr); % Concatenating
    Phrelngr = horzcat(Phrelngr,P_rel_ngr); % Concatenating

    % FORMING THE ABSOLUTE CONSTRAINTS of NORTH GRADIENTS
    coef_abs_ngr = opt.stat(istat).coef_abs_ngr;
    
        
        for inter=1:(t_ngr+na)
            Habsngr(inter, inter) = +1;  % design matrix for the east gradients pseudo observtaions as constraints
            Phabsngr(inter,inter) = 1./(coef_abs_ngr)^2; % [1/cm^2] weight matrix of the design matrix for the North gradients constraints
        end
        
    
    clear H_rel_ngr P_rel_ngr 
end

% FORMING THE O-C VECTOR FOR THE RELATIVE & ABSOLUTE CONSTRAINTS
och_rel(size(Hrelngr,1),1) = 0;  % O-C vector for the relative constraints of  North gradients
och_abs(size(Habsngr,1),1) = 0;  % O-C vector for the absolute constraints of  North gradients


%..........................................................................

% ELIMINATE PARAMETERS ACCORDING TO opt

for istat=na:-1:1

   if opt.stat(istat).ngr_inc == 0
        
        if opt.constr_rel_ngr == 1 && opt.constr_abs_ngr == 0
            Hrelngr(:,sum(num_inter_ngr(1:istat-1))+istat:sum(num_inter_ngr(1:istat))+istat) = [];
            Hrelngr(sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat)),:) = [];
            
            Phrelngr(:,sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat))) = [];
            Phrelngr(sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat)),:) = [];
            
            och_rel(sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat))) = [];
           
        end

        if opt.constr_rel_ngr == 0 && opt.constr_abs_ngr == 1
            Habsngr(:,sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat))) = [];
            Habsngr(sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat)),:) = [];
            
            Phabsngr(:,sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat))) = [];
            Phabsngr(sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat)),:) = [];
            
            och_abs(sum(num_inter_ngr(1:istat-1))+istat:sum(num_inter_ngr(1:istat))+istat) = [];
        end
        
        if opt.constr_rel_ngr == 1 && opt.constr_abs_ngr == 1            
            Hrelngr(:,sum(num_inter_ngr(1:istat-1))+istat:sum(num_inter_ngr(1:istat))+istat) = [];
            Hrelngr(sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat)),:) = [];
           
            Phrelngr(:,sum(num_inter_ngr(1:istat-1)):sum(num_inter_ngr(1:istat))) = [];
            Phrelngr(sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat)),:) = [];
           
            och_rel(sum(num_inter_ngr(1:istat-1)):sum(num_inter_ngr(1:istat))) = [];

            Habsngr(:,sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat))) = [];
            Habsngr(sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat)),:) = [];
            
            Phabsngr(:,sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat))) = [];
            Phabsngr(sum(num_inter_ngr(1:istat-1))+1:sum(num_inter_ngr(1:istat)),:) = [];
            
            och_abs(sum(num_inter_ngr(1:istat-1))+istat:sum(num_inter_ngr(1:istat))+istat) = [];
        end
       
    end



end




%..........................................................................
% TAKE RELATIVE OR ABSOLUTE ACCORDING TO opt
if opt.constr_rel_ngr==1 && opt.constr_abs_ngr==0
   
    H4=Hrelngr;
    Ph4=Phrelngr;
    och4=och_rel;
   
end

if opt.constr_rel_ngr==0 && opt.constr_abs_ngr==1
  
    H4=Habsngr;
    Ph4=Phabsngr;
    och4=och_abs;
   
end

if opt.constr_rel_ngr==1 && opt.constr_abs_ngr==1
   
    H4=vertcat(Hrelngr,Habsngr);
    Ph4=blkdiag(Phrelngr,Phabsngr);
    och4=vertcat(och_rel,och_abs);
   
end

if opt.constr_rel_ngr==0 && opt.constr_abs_ngr==0
   
    if col_A4==0
       H4=[];
    else H4(1,col_A4)=0;
    end    
    Ph4=[];
    och4=[];
   
end

