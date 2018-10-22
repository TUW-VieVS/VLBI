% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   constraints of antenas' coordinates
%
%   Reference: 
%
%   Input:	
%       'opt'                  structure array         (for info. /DOC/opt.doc)
%       'na'                   (1,1)                   number of antennas
%
%   Output:
%       'HX'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'PhX'         structure array     weight matrix of the pseudo-observation equations 
%       'ochX'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   July 2012 by Claudia Tierno Ros
%
%   Revision: 
%  
% ************************************************************************


function [H13,H14,H15,Ph13,Ph14,Ph15,och13,och14,och15] = hpoc_xyz_scan(opt,na)

if opt.pw_stc==0
   
   % X
    H13(na-1,na) = 0;
    Ph13(na-1,na-1) = 0;
    och13=zeros(size( Ph13,1),1);
   
   % Y
    H14=H13;
    Ph14=Ph13;
    och14=zeros(size( Ph14,1),1);
    
   % Z
    H15=H13;
    Ph15=Ph13;
    och15=zeros(size( Ph15,1),1);
    
else if opt.pw_stc==1

   % X         
    H13=zeros(na-1,na);
    for k=1:length(na)
        H13(k,k)=1;
        H13(k,k+1)=-1;
        for j=1:na
            Ph13(k,k)=1./opt.stat(k).coef_xyz^2; % [1/cm^2] weight matrix of the design matrix for the coordinates constraints
        end
     end
    och13=zeros(size( Ph13,1),1);
 
   % Y
    H14=H13;
    Ph14=Ph13;
    och14=och13;
  
   % Z
    H15=H13;
    Ph15=Ph13;
    och15=och13;
  
    end
end


if opt.constr_xyz == 0
    
   % X
    H13(na-1,na) = 0;
    Ph13(na-1,na-1) = 0;
    och13(size( Ph13,1),1)=[];
    
   % Y
    H14=H13;
    Ph14=Ph13;
    och14(size( Ph14,1),1)=[];
    
   % Z
    H15=H13;
    Ph15=Ph13;
    och15(size( Ph15,1),1)=[];
end

end


