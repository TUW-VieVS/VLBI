% ************************************************************************
%   Description:
%       This function creates the B matrix - NNT/NNR condition on
%       coordinates/velocities
%
%
%
%   Input:
%      x0           reduced a priori coordinates
%      y0
%      z0
%      ln           number of estimated stations
%      excidant     indices of excluded stations from the NNT/NNR condition
%      
%
%
%   Output: 
%      B            matrix with coefficients for free-network constraint
%
%
%   External calls: 	
%      -               					    											
%       
%   Coded for VieVS: 
%   03 Aug 2010 by Hana Spicakova
%
%   Revision: 
%   11 Oct 2010 by Hana Spicakova: scale added (now it can be choosen between 6- and 7-parameter transformation)
%
%%



function B=nntnnr(x0,y0,z0,ln,excidant,scale)
   
    clear B
    
    if scale=='0'
        for j=1:ln
            Bj =[  1      0      0 
                   0      1      0
                   0      0      1
                   0    -z0(j)  y0(j)
                  z0(j)   0    -x0(j)
                 -y0(j)  x0(j)   0];

            B(:,j)      =Bj(:,1); % X
            B(:,j+ln)   =Bj(:,2); % Y
            B(:,j+ln*2) =Bj(:,3); % Z
        end
        
    elseif scale=='1'
        for j=1:ln
            Bj =[  1      0      0 
                   0      1      0
                   0      0      1
                   0    -z0(j)  y0(j)
                  z0(j)   0    -x0(j)
                 -y0(j)  x0(j)   0
                  x0(j)  y0(j)  z0(j)];

            B(:,j)      =Bj(:,1); % X
            B(:,j+ln)   =Bj(:,2); % Y
            B(:,j+ln*2) =Bj(:,3); % Z
        end
    end
    
   
    if isempty(excidant)==0
        if scale=='0'
            k=6;
        elseif scale=='1'
            k=7;
        end
        
        for i=1:length(excidant)
            B(:,excidant(i))          =zeros(1,k);
            B(:,excidant(i)+ln)       =zeros(1,k);
            B(:,excidant(i)+ln*2)     =zeros(1,k);
        end
    end



