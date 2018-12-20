% ************************************************************************
%   Description:
%   Function to create the N and b matrixes from H, Ph and och
%
%   Reference: 
%
%   Input:	
%       'HX'          structure array      design matrix of the pseudo-observation equations (constraints) of constraint 'x'
%       'PhX'         structure array      weight matrix of the pseudo-observation equations of constraint 'x'
%       'ochX'        structure array      o-c vector of constraints (zero vector) of constraint 'x'
%       'Nh'          matrix               N matrix of constraints
%       'bh'          matrix               b matrix of constraints 
%       'ncol'        matrix               num. of column for the next matrix to start
%       'col_Ablk'    (1,1)                number of columns in Ablk
%       'Ngc'         matrix     N matrix of constraints for sinex output
%       'bgc'         matrix     b matrix of constraints for sinex output
%       'Hblk'        matrix     H matrix of all the constraints
%       'och_total'   matrix     o-c vector of all the constraints
%
%   Output:
%       'Nh'          matrix     N matrix of constraints
%       'bh'          matrix     b matrix of constraints 
%       'ncol'        matrix     num. of column for the next matrix to start
%       'Ngc'         matrix     N matrix of constraints for sinex output
%       'bgc'         matrix     b matrix of constraints for sinex output
%       'Hblk'        matrix     H matrix of all the constraints
%       'och_total'   matrix     o-c vector of all the constraints
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   27 Aug 2012 by Claudia Tierno Ros
%
%   Revision: 
%
% ************************************************************************


 function [Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(Hx,Phx,ochx,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total)
 

temp=size(Hx,2);
Hx(~any(Hx,2),:) = [];

 if ~isempty(Hx)%isempty(Hx)==0
   
   % Making H the correct size  
   H=zeros(size(Hx,1),col_Ablk);
   H(1:size(Hx,1),ncol+1:ncol+size(Hx,2))=Hx;

   N=sparse(H'*Phx*H);   % N of constraint 'x'                
 
   Hblk=sparse(vertcat(Hblk,H));              %Needed for calculating the residuals
   och_total=sparse(vertcat(och_total,ochx)); %Needed for calculating the residuals

   if ~isempty(ochx) % b of constraint 'x'
      b=H'*Phx*ochx; 
   end

   Nh=Nh+N;  % N of all the constraints

   if exist('b','var')  % b of all the constraints
    bh=bh+b;
   end
   
   
   % Global solution & sinex output
    if  opt.ascii_snx == 1 || opt.global_solve==1
    
      Hglob=sparse(horzcat(H,zeros(size(H,1),col_A_add)));
%      Hglob(:,col_est)=0; % remove constraints of parameters which will be written into SINEX file
      Nglob_const=sparse(Hglob'*Phx*Hglob);
    
      if ~isempty(ochx)
        bglob_const=Hglob'*Phx*ochx; 
      end
     
      Ngc=sparse(Ngc+Nglob_const);
     
      if exist('bglob_const','var') 
         bgc=bgc+bglob_const;
      end
     
    end
 end
ncol=ncol+temp;
    
    

