% ************************************************************************
%   Description:
%   function to delete reference clock 
%
%   Reference: 
%
%   Input:	
%       'A_session'       structure array        design matrix of the real-observation equations
%       'num_inter_clk'   (1, num. sources)      number of estimation intervals for clocks
%       'nistat'          (1,1)                  number of reference station
%
%   Output:
%       'A_session'     structure array     design matrix of the real-observation equations
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   27 Aug 2012 by Claudia Tierno Ros
%
%   Revision: 
%
% ************************************************************************





function [A_session] = delref_scan_A(A_session,num_inter_clk,nistat,opt)

    % Eliminate corresponding columns
    
    if nistat==1 % if reference station is station number 1
        A_session(1).sm(:,1:1+num_inter_clk(1)) = []; % from piecewise linear functions
        if opt.pw_clk~=2
            A_session(2).sm(:,1:2)=[]; % from quadratic clock functions
        else
            A_session(2).sm(:,1)=[];
        end
    else
        % if reference station is diferent than number 1
        A_session(1).sm(:,sum(num_inter_clk(1:nistat-1))+nistat:(sum(num_inter_clk(1:nistat-1))+nistat-1)+num_inter_clk(nistat)+1) = []; % from piecewise linear functions
       if opt.pw_clk~=2
        A_session(2).sm(:,nistat*2-1:nistat*2) = []; % from quadratic clock functions
       else
        A_session(2).sm(:,nistat)=[];
       end
  
    
    end
    
    



