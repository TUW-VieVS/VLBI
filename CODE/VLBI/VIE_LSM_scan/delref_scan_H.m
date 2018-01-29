% ************************************************************************
%   Description:
%   function to delete reference clock from design matrix, weight matrix,
%   and o-c vector of constraints
%
%   Reference: 
%
%   Input:	
%       'HX'               structure array      design matrix of the pseudo-observation equations (constraints)
%       'PhX'              structure array      weight matrix of the pseudo-observation equations 
%       'ochX'             structure array      o-c vector of constraints (zero vector)
%       'opt'              structure array      (for info. /DOC/opt.doc)
%       'num_inter_clk'   (1, num. sources)     number of estimation intervals for clocks
%       'nistat'          (1,1)                 number of reference station
%
%   Output:
%       'HX'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'PhX'         structure array     weight matrix of the pseudo-observation equations 
%       'ochX'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   27 Aug 2012 by Claudia Tierno Ros
%
%   Revision: 
%
% ************************************************************************

function [H1,Ph1,och1,H2,Ph2] = delref_scan_H(opt,H1,Ph1,H2,Ph2,och1,num_inter_clk,nistat,na,num_inter_qclk)

    % Eliminate corresponding columns
    
    if nistat==1
       
        H1(:,1:1+num_inter_clk(1)) = [];
        H1(1:num_inter_clk(1),:)=[];
        
        if opt.pw_clk~=2
        H2(:,1:2)=[];
         H2(1,:)=[];
        else
        H2(:,1)=[];
        end
        Ph2(:,nistat)=[];
        Ph2(nistat,:)=[];
       
        
        Ph1(:,1:num_inter_clk(1))=[];
        Ph1(1:num_inter_clk(1),:)=[];
        
                
        if opt.constr_clk==1
        och1(1:num_inter_clk(1),:)= [];
        end
        
    else

    H1(:,sum(num_inter_clk(1:nistat-1))+nistat:(sum(num_inter_clk(1:nistat-1))+nistat-1)+num_inter_clk(nistat)+1)=[];
    H1(sum(num_inter_clk(1:nistat-1))+1:(sum(num_inter_clk(1:nistat-1))-1)+num_inter_clk(nistat)+1,:)=[];

    
    
    sum_qclk(1:na+1)=0;
    sum_qclk(2:na+1)=cumsum(num_inter_qclk);
    vecm = 0:1:na;
    summ=sum_qclk-vecm;
    
    if opt.pw_clk~=2
    H2(:,nistat*2-1:nistat*2)=[];
        H2(nistat,:)=[];
    Ph2(:,nistat)=[];
    Ph2(nistat,:)=[];
    else
    H2(:,nistat)=[];
    H2(summ(nistat)+1:summ(nistat+1),:) = [];
    Ph2(summ(nistat)+1:summ(nistat+1),:) = [];
    Ph2(:,summ(nistat)+1:summ(nistat+1)) = [];

    end
    
    Ph1(:,sum(num_inter_clk(1:nistat-1))+1:(sum(num_inter_clk(1:nistat-1))-1)+num_inter_clk(nistat)+1)=[];
    Ph1(sum(num_inter_clk(1:nistat-1))+1:(sum(num_inter_clk(1:nistat-1))-1)+num_inter_clk(nistat)+1,:)=[];
    
    
    if opt.constr_clk==1
        
        och1(sum(num_inter_clk(1:nistat-1))+1:(sum(num_inter_clk(1:nistat-1))-1)+num_inter_clk(nistat)+1,:)=[];
        
    end
    end
    



