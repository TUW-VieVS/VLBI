% ************************************************************************
%   Description:
%   function to form NNT/NNR condition equations for datum definition of
%   TRF
%
%   Reference: 
%   nnt_cond.m
%
%   Input:	
%       'na'         (1,1)                      number of antennas
%       'x0'         (1,na)                     apriori TRF (x) coordinates of all antennas in the session
%       'y0'         (1,na)                     apriori TRF (y) coordinates of all antennas in the session
%       'z0'         (1,na)                     apriori TRF (z) coordinates of all antennas in the session
%       'opt'        structure array            (for info. /DOC/opt.doc)
%       'sum_dj'     (1,number of models)       total number of estimates for each included model
%       'N'          (sum_dj(end),sum_dj(end))  datum free normal equation matrix         
%
%   Output:
%       'N' (sum_dj(end)+cond,sum_dj(end)+cond) normal equation matrix with NNT/NNR condition equations        
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   31 July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************

function [N] = nnt_cond_scan(na,x0,y0,z0,opt,sum_dj,N) 
             
        Bx = []; By = []; Bz = []; B2 = [];
    
        xs = mean(x0); ys = mean(y0); zs = mean(z0);
        xi = x0 - xs; yi = y0 - ys; zi = z0 - zs;
        cc = 1/sqrt(x0*x0'+y0*y0'+z0*z0');
        xii = cc*x0; yii = cc*y0; zii = cc*z0;
        n_xyz=1; % number of estimates 
                 % Always 1 because coord. always estimated from the whole
                 % session 
                 
        for istat = 1 : na  
    
        
           
        %--------------------
        btx(3,1) = 0; bty = btx; btz = btx;
        bsx(1,1) = 0; bsy = bsx; bsz = bsx;
        brx(3,1) = 0; bry = brx; brz = brx;
        %--------------------
        btx(1,1) = ones(1,1); 
        bsx(1,1) = xii(istat)*ones(1,1); 
        brx(2,1) = zii(istat)*ones(1,1); 
        brx(3,1) =-yii(istat)*ones(1,1); 
        %--------------------
        bty(2,1) = ones(1,1); 
        bsy(1,1) = yii(istat)*ones(1,1);
        bry(1,1) =-zii(istat)*ones(1,1); 
        bry(3,1) = xii(istat)*ones(1,1); 
        %--------------------
        btz(3,1) = ones(1,1); 
        bsz(1,1) = zii(istat)*ones(1,1);
        brz(1,1) = yii(istat)*ones(1,1); 
        brz(2,1) =-xii(istat)*ones(1,1); 
        %--------------------
        if opt.stat(istat).nnt_inc == 0;
           clear btx, btx(3,1) = 0;
           clear bty, bty = btx;
           clear btz, btz = btx; 
        end
        if opt.stat(istat).nnr_inc == 0;
           clear brx, brx(3,1) = 0;
           clear bry, bry = brx; 
           clear brz, brz = brx; 
        end
        if opt.stat(istat).nns_inc == 0;
           clear bsx, bsx(1,1) = 0; 
           clear bsy, bsy = bsx; 
           clear bsz, bsz = bsx; 
        end
        
        bx = vertcat(btx,bsx,brx);
        by = vertcat(bty,bsy,bry);
        bz = vertcat(btz,bsz,brz);
        
        Bx = horzcat(Bx,bx); By = horzcat(By,by); Bz = horzcat(Bz,bz);
        
        clear btx bty btz brx bry brz bsx bsy bsz bx by bz
       end    
    
      
   B1 = horzcat(Bx,By,Bz); 
   B1(~any(B1,2),:) = [];
   B2(size(B1,1),sum_dj(13)) =  0; 
   B = horzcat(B2,B1); 
   K(size(B,1),size(B,1)) = 0; 
   N2 = horzcat(vertcat(N,B),vertcat(B',K)); 
   clear N
   N = N2; 
      
    