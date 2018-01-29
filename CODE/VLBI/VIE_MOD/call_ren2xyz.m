% ************************************************************************
%   Description:
%   Transformation into geocentric Cartesian system. 
% 
%   Input:										
%      matm                array with time and correction data
%                          (mjd r e n) in [m]
%      ant                 station cartesian coordinates from catalogue 
%                
%   Output:
%      mjdxyz              mjd and xyz corrections [m]
%
%   Coded for VieVS: 
%   27 March 2012 by Hana Spicakova
%
%   Revision: 
% *************************************************************************
function [mjdxyz]=call_ren2xyz(matm,ant)

   % geodetic latitude
    phi1=cart2phigd(ant); 
    lam1=atan2(ant(2),ant(1));

    ldr=size(matm,1);
    phi(1:ldr,1)=phi1;
    lam(1:ldr,1)=lam1;

    % Correction to station vector
    dren=matm(:,2:4);
    [dxyz]=ren2xyz(dren,phi,lam);    % [m]
    
    mjdxyz=[matm(:,1) dxyz];
    


   
   
   
   
   
   
   
   
   