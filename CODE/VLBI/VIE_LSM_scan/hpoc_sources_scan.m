% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   constraints of sources' coordinates
%
%   Reference: 
%
%   Input:	
%       'num_inter_sou'        (1,1)                   number of estimation intervals source coord.
%       'opt'                  structure array         (for info. /DOC/opt.doc)
%       'ns'                   (1,1)                   number of sources
%
%   Output:
%       'HX'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'PhX'         structure array     weight matrix of the pseudo-observation equations 
%       'ochX'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   29 August 2012 by Claudia Tierno Ros
%
%   Revision: 
%   18 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with
%   NNR condition
%   17 Jul 2014 by Hana Krasna: very loose constrains of 10 mas are applied
%   on all sources together with the NNR condition to prevent the matrix
%   being singular
% ************************************************************************

function [H11,H12,Ph11,Ph12,och11,och12] = hpoc_sources_scan(ns,opt,num_inter_sou)

Hra = []; Phra = []; Hde = []; Phde = []; 
t_sou=sum(num_inter_sou); %total estimate interval for sources


for isou = 1 : ns 
    % FORMING THE CONSTRAINTS of SOURCES COORDINATES  
    if opt.pw_sou == 1 % piecewise

        coef_rade = opt.source(isou).coef_rade;
        int_rade = opt.source(isou).int_rade;
        H_ra(t_sou,num_inter_sou(isou)+1) = 0; 
        P_hra(t_sou,num_inter_sou(isou))= 0;        

        sumrade = sum(num_inter_sou(1:isou-1));

        for inter = 1:num_inter_sou(isou) 
            H_ra(sumrade + inter, inter) = +1;  % design matrix for the right ascension pseudo observtaions as constarints
            H_ra(sumrade + inter, inter+1) = -1; % design matrix for the declination pseudo observtaions as constarints
            P_hra(sumrade + inter, inter) = 1./coef_rade^2; % [1/mas^2] weight matrix of the design matrix for the sources coordinates constraints
        end     
        Hra = horzcat(Hra,H_ra); % Concatenating
        Phra = horzcat(Phra,P_hra); % Concatenating
        Hde = Hra; 
        Phde = Phra; 
        clear H_ra P_hra  
    end
end

p=1/10^2; % constraints 10 mas
if opt.est_sourceNNR==1
    Hra=diag(ones(ns,1));
    Hde=diag(ones(ns,1));
    Phra=eye(ns).*p;
    Phde=eye(ns).*p;
end



% % -------------------------------------------------------------------------
% % FORMING THE O-C VECTOR FOR THE CONSTRAINTS
% oc_hra(size(Hra,1),1) = 0;  % O-C vector for the constraints of  right ascension coordinates
% oc_hde(size(Hde,1),1) = 0;  % O-C vector for the constraints of  declination coordinates


%--------------------------------------------------------------------------
% ELIMINATE SOURCES ACCORDING TO opt
if opt.pw_sou == 1 & opt.constr_sou == 0
    for isou=ns:-1:1
        if   opt.source(isou).rade_inc == 0 
             if isou==1
                Hra(:,1:num_inter_sou(1)+1) =[];
                Hra(1:num_inter_sou(1),:) =[];

                Hde(:,1:num_inter_sou(1)+1) =[];
                Hde(1:num_inter_sou(1),:) =[];

                Phra(:,1:num_inter_sou(1)) =[];
                Phra(1:num_inter_sou(1),:) =[];

                Phde(:,1:num_inter_sou(1)) =[];
                Phde(1:num_inter_sou(1),:) =[];

             else
                Hra(:,sum(num_inter_sou(1:isou-1))+isou:sum(num_inter_sou(1:isou))+isou) = [];
                Hra(sum(num_inter_sou(1:isou-1))+1:sum(num_inter_sou(1:isou)),:) = [];

                Hde(:,sum(num_inter_sou(1:isou-1))+isou:sum(num_inter_sou(1:isou))+isou) = [];
                Hde(sum(num_inter_sou(1:isou-1))+1:sum(num_inter_sou(1:isou)),:) = [];

                Phra(:,sum(num_inter_sou(1:isou-1))+1:sum(num_inter_sou(1:isou))) = [];
                Phra(sum(num_inter_sou(1:isou-1))+1:sum(num_inter_sou(1:isou)),:) = [];

                Phde(:,sum(num_inter_sou(1:isou-1))+1:sum(num_inter_sou(1:isou))) = [];
                Phde(sum(num_inter_sou(1:isou-1))+1:sum(num_inter_sou(1:isou)),:) = [];
             end

        end
    

%        if isou==1
%            oc_hra(1:num_inter_sou(1))=[];
%        else
%            oc_hra(sum(num_inter_sou(1:isou-1))+1:sum(num_inter_sou(1:isou))) = [];
%            %oc_hde = oc_hra; 
%        end
%   
%     oc_hde = oc_hra; 
    end
end
    
    

och11(1:size(Hra,1),1)=0;
och12=och11;

H11 = Hra; H12 = Hde; 
Ph11 = Phra; Ph12 = Phde; 
%och11 = oc_hra; 
%och12 = oc_hde;
