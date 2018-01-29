% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   local ties (station coordinates) at the co-location sites
%
%   Reference: 
%
%   Input:	
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%       'n_'         structure array     number of estimates (pwlo or one offset)
%       'opt'        structure array     (for info. /DOC/opt.doc)
%       'na'         (1,1)               number of antennas
%       't_'         structure array     total estimate intervals of a specific estimate in a session
%       'antenna'    structure array     antenna information of the session
%
%   Output:
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%       't_'         structure array     total estimate intervals of a specific estimate in a session
%       'nlt'        (1,1)               number of local ties
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   08 Jun 2015 by Younghee Kwak
%
%   Revision: 
%   17 Feb 2016 by Younghee Kwak modified for the standard VieVS version
% ************************************************************************
function [H,Ph,och,t_,nlt] = hpoc_lt(H,Ph,och,n_,opt,na,t_,antenna)

%-----------------------------------------%
% should be a separate module or added to the superstation file %
% now only hobart antennas % There might be twin/sibling telescopes in Wettzell, Onsala, Hartebeesthoek, Yebes, Kashima ...
%fid=fopen('local_tie.dat','r');
%fgetl(fid);  % skip the first line
%ista=0;
%while ~feof(fid)
%    str=fgetl(fid);
%    ista=ista+1;
%    delta_X(ista)  =str2num(str(22:33));  sigma_X(ista)  =str2num(str(34:42));  % m
%    delta_Y(ista)  =str2num(str(43:54));  sigma_Y(ista)  =str2num(str(55:63));  % m
%    delta_Z(ista)  =str2num(str(64:75));  sigma_Z(ista)  =str2num(str(76:84));  % m
%end
%fclose(fid);
nlt=1; % number of local ties % It should be counted from local tie information
for istat = 1:na
    if antenna(istat).name == 'HOBART26'
        ncombi(1,1)=istat;
    end
    if antenna(istat).name == 'HOBART12'
        ncombi(1,2)=istat;
    end
%    if antenna(istat).name == 'HARTRAO '
%        ncombi(2,1)=istat;
%    end
%    if antenna(istat).name == 'HART15M '
%        ncombi(2,2)=istat;
%    end
end 
delta_X(1)  =  246.5743 ;  % local tie between two reference points % m
delta_Y(1)  =   73.5364 ;  
delta_Z(1)  = -146.1538 ;  
%-----------------------------------------%
sessionTimeConsumption=max([antenna.lastObsMjd])-min([antenna.firstObsMjd]);    % needed later when writing epochs where parameters are valid 
if sessionTimeConsumption>0.9
    midmjd=floor(max([antenna.lastObsMjd])-0.5)+0.5;        %hana
else %session lasts less than one (0.9) day, eg. one hour
    midmjd=floor(max([antenna.lastObsMjd]))+0.5;          % hana jun11
end

Hcomx = []; Phcomx = []; Hcomy = []; Phcomy = []; Hcomz = []; Phcomz = [];
ochcomx = []; ochcomy = []; ochcomz = [];
t_.xyz = sum([n_.xyz]); % total estimate interval of coordinates in a session
ncolo=1; % number of co-location sites % currently only 1, hobart
for icom = 1:ncolo % number of co-location site
    sumxyzj = 0;
    sumxyzk = 0;

    j=ncombi(icom,1);  % istat for co-location site  !! the order at antenna.mat
    k=ncombi(icom,2);  % istat for co-location site  !! the order at antenna.mat

    antenna(j).aprX=antenna(j).x+antenna(j).vx*(midmjd-antenna(j).epoch)/365.25;
    antenna(j).aprY=antenna(j).y+antenna(j).vy*(midmjd-antenna(j).epoch)/365.25;
    antenna(j).aprZ=antenna(j).z+antenna(j).vz*(midmjd-antenna(j).epoch)/365.25;
    antenna(k).aprX=antenna(k).x+antenna(k).vx*(midmjd-antenna(k).epoch)/365.25;
    antenna(k).aprY=antenna(k).y+antenna(k).vy*(midmjd-antenna(k).epoch)/365.25;
    antenna(k).aprZ=antenna(k).z+antenna(k).vz*(midmjd-antenna(k).epoch)/365.25;


    H_comx(icom).h(n_(j).xyz,t_.xyz) = 0; H_comy = H_comx; H_comz = H_comx;        % initialization for constraint matrix        !! the matrix size: (1) row    = no. of co-location sites * number of x per
                                                                                   %                                                                 (2) column = no. of x parameters 
    Ph_comx(icom).h(n_(j).xyz,n_(j).xyz)= 0; Ph_comy = Ph_comx; Ph_comz = Ph_comx; % initialization for constraint weight matrix !! diagonal square matrix with same number of row of H_comx               
    oc_hcomx(icom).h(n_(j).xyz,1) = 0; oc_hcomy(icom).h(n_(j).xyz,1) = 0; oc_hcomz(icom).h(n_(j).xyz,1) = 0;


    for i = 1:j-1 % for both antennas j & k
        sumxyzj = sumxyzj + n_(i).xyz;  % counting number of x estimates right before the estimate to use
                                        % if it is the 6th site, count x numbers till end of 5th site
    end
    for i = 1:k-1 
        sumxyzk = sumxyzk + n_(i).xyz;
    end

    for inter = 1:n_(j).xyz       % for convenience, assume the numbers of x are same for both antennas n_(j).xyz = n_(k).xyz
        H_comx(icom).h(inter,sumxyzj+inter) = -1;  %  e.g. x of 'HOBART26(Hb)' % 1*X_Hb - 1*X_Ho = Delta X - (X_apr_Ho - X_apr_Hb) 
        H_comx(icom).h(inter,sumxyzk+inter) = +1;  %  e.g. x of 'HOBART12(Ho)'
        Ph_comx(icom).h(inter,inter) = 1./...
                        (3^2); %  [1/cm^2] 3cm weight matrix of the design matrix for the local tie
        oc_hcomx(icom).h(inter,1) = (delta_X(icom)-(antenna(k).aprX-antenna(j).aprX))*100;  % local tie vector  m -> cm

        H_comy(icom).h(inter,sumxyzj+inter) = -1;
        H_comy(icom).h(inter,sumxyzk+inter) = +1;
        Ph_comy(icom).h(inter,inter) = 1./...
                        (3^2);
        oc_hcomy(icom).h(inter,1) = (delta_Y(icom)-(antenna(k).aprY-antenna(j).aprY))*100;

        H_comz(icom).h(inter,sumxyzj+inter) = -1; 
        H_comz(icom).h(inter,sumxyzk+inter) = +1; 
        Ph_comz(icom).h(inter,inter) = 1./...
                        (3^2);
        oc_hcomz(icom).h(inter,1) = (delta_Z(icom)-(antenna(k).aprZ-antenna(j).aprZ))*100;
    end

    % Concatenating
    Hcomx = vertcat(Hcomx,H_comx(icom).h);
    Phcomx = blkdiag(Phcomx,Ph_comx(icom).h);
    ochcomx = vertcat(ochcomx,oc_hcomx(icom).h);
    Hcomy = vertcat(Hcomy,H_comy(icom).h);
    Phcomy = blkdiag(Phcomy,Ph_comy(icom).h);
    ochcomy = vertcat(ochcomy,oc_hcomy(icom).h);
    Hcomz = vertcat(Hcomz,H_comz(icom).h);
    Phcomz = blkdiag(Phcomz,Ph_comz(icom).h);
    ochcomz = vertcat(ochcomz,oc_hcomz(icom).h);
end

%save h.mat H_comx;   % test
%save ph.mat Ph_comx; %test

% add to original constrain matrix
H(13).sm   = vertcat(H(13).sm,   Hcomx);   
Ph(13).sm  = blkdiag(Ph(13).sm,  Phcomx); 
och(13).sv = vertcat(och(13).sv, ochcomx); 
H(14).sm   = vertcat(H(14).sm,   Hcomy);   
Ph(14).sm  = blkdiag(Ph(14).sm,  Phcomy); 
och(14).sv = vertcat(och(14).sv, ochcomy); 
H(15).sm   = vertcat(H(15).sm,   Hcomz);   
Ph(15).sm  = blkdiag(Ph(15).sm,  Phcomz); 
och(15).sv = vertcat(och(15).sv, ochcomz); 
