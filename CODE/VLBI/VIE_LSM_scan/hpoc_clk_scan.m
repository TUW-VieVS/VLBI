% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   constraints of clocks
%
%   Reference: 
%
%   Input:	
%       'num_inter_clk'        (1,num. of antennas)     number of estimation intervals (pwlo or one offset)
%       'opt'                  structure array         (for info. /DOC/opt.doc)
%       'na'                   (1,1)                    number of antennas
%       'c'                    (1,1)                    velocity of light in vacuum  
%       'col_A2'               (1,1)                    number of columns of A_session(2)
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

function [H1,Ph1,H2,Ph2,och1] = hpoc_clk_scan(num_inter_clk,opt,na,c,col_A2)

H1 = []; Ph1 = [];
t_clk = sum(num_inter_clk); % total estimate interval of CLOCKs in a session
for istat = 1:na
    % FORMING THE CONSTRAINTS of CLOCKSs
    coef_clk = opt.stat(istat).coef_clk; 
    int_clk = opt.stat(istat).int_clk;
     H_clk(t_clk,num_inter_clk(istat)+1) = 0;
    Ph_clk(t_clk,num_inter_clk(istat))= 0;
        sumclk = 0;
        for i = 1:istat-1
            sumclk = sumclk + num_inter_clk(i);
        end
        for inter = 1:num_inter_clk(istat) 
            H_clk(sumclk+inter,inter) = +1;  % design matrix for the clock pseudo observations as constarints
            H_clk(sumclk+inter,inter+1) = -1; % design matrix for the clock pseudo observations as constarints
            Ph_clk(sumclk+inter,inter) = 1./coef_clk^2; % [1/cm^2] weight matrix of the design matrix for the clocks constraints
             
        end                     
    H1 = horzcat(H1,H_clk); % Concatenating
    Ph1 = horzcat(Ph1,Ph_clk); % Concatenating
    clear H_clk Ph_clk
end
% FORMING THE O-C VECTOR FOR THE CONSTRAINTS
och1(size(H1,1),1) = 0;  % O-C vector for the clock constraints

if opt.constr_clk == 0
    H1(1:size(H1,1),1:size(H1,2)) = 0; 
    Ph1(1:size(Ph1,1),1:size(Ph1,2)) = 0;
    och1 = [];
end


H2(col_A2-na+1,col_A2) = 0;

Ph2(col_A2-na+1,col_A2-na+1) = 0;

end

