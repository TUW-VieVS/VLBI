% ************************************************************************
%   Description:
%       This function creates the B matrix - NNR condition on source
%       coordinates
%
%
%
%   Input:
%      RA           a priori coordinates of sources [rad,rad]
%      De
%      lns          number of estimated sources
%      excidsouc    indices of excluded sources from the NNR condition
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
%   04 Oct 2010 by Hana Spicakova:  choose between NNR+dz (3 rotations and
%                                   translation in declination) and only NNR
%%

function B=nnr_source(RA,De,lns,excidsouc,sou_dz)
    clear B

    if sou_dz=='1'
        for j=1:lns
            Bj =[ tan(De(j))*cos(RA(j))  -sin(RA(j))
                  tan(De(j))*sin(RA(j))   cos(RA(j)) 
                  -1                      0
                   0                      1];

            B(:,j)       =Bj(:,1); % RA
            B(:,j+lns)   =Bj(:,2); % De
        end
      
    elseif sou_dz=='0'
        for j=1:lns
            Bj =[ tan(De(j))*cos(RA(j))  -sin(RA(j))
                  tan(De(j))*sin(RA(j))   cos(RA(j)) 
                  -1                      0];

            B(:,j)       =Bj(:,1); % RA
            B(:,j+lns)   =Bj(:,2); % De
        end
    end

    if isempty(excidsouc)==0
        if sou_dz=='1'
            k=4;
        elseif sou_dz=='0'
            k=3;
        end
        for i=1:length(excidsouc)
            B(:,excidsouc(i))          =zeros(1,k);
            B(:,excidsouc(i)+lns)      =zeros(1,k);
        end
    end
    
