% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   constraints of common zenith wet delays at the co-location sites
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
%       'c'          (1,1)               velocity of light in vacuum  
%       'antenna'    structure array     antenna information of the session
%       't'          structure array     estimation intervals of clk, zwd, ngr, egr, xyz
%
%   Output:
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%       't_'         structure array     total estimate intervals of a specific estimate in a session
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   based on Kamil's hpoc_* codes
%   14 May 2015 by Younghee Kwak
%   
%   Revision: 
%   17 Feb 2016 by Younghee Kwak modified for the standard VieVS version
%   01 Feb 2017 by Younghee Kwak added 3 more sites: Hobart + Hartebeesthoek, Wettzell, Yebes
%                                check co-located sites
%   12 Feb 2017 by Younghee Kwak read proper trop. tie data
%   13 Feb 2017 by Younghee Kwak skip the non-matching intervals 
% ************************************************************************
function [H,Ph,och,t_] = hpoc_comzwd(H,Ph,och,n_,opt,na,c,antenna,t)

%-----------------------------------------%
% should be a separate module or added to the superstation file %
% four co-located VLBI sites with twin/sibling telescopes, i.e. Hobart, Hartebeesthoek, Wettzell, Yebes
fid=fopen('tropo_tie.dat','r');
fgetl(fid);  % skip the first line
ista=0;
while ~feof(fid)
    str=fgetl(fid);
    ista=ista+1;
    station_name_new{ista}=str(1:8);    % 8 characters vgos antenna
    station_name_old{ista}=str(10:17);  % 8 characters legacy antenna
    %height_diff(ista)=str2num(str(12:19));   % m
    delta_ZTD_f(ista)  =str2num(str(53:60));   % mm
    delta_ZHD_f(ista)  =str2num(str(65:72));   % mm
    delta_ZWD_f(ista)  =str2num(str(76:83))*0.1;   % mm --> cm
    nsta=ista;
end
fclose(fid);
%delta_ZWD(1)  = 1.0 ; % tropo-tie (ZWD) between two reference points  [cm] % ZWD_Hb - ZWD_Ho ~ 1cm

icolo=0;
for istat = 1:na
    if antenna(istat).name == 'HOBART26'
        for jstat = 1:na
            if antenna(jstat).name == 'HOBART12'
                % if both antennas are included in the list, treat them as a co-located site
                icolo=icolo+1;
                ncombi(icolo,1)=istat;
                ncombi(icolo,2)=jstat;
                for ista=1:nsta 
                    if (strcmp(antenna(jstat).name, station_name_new(ista)) && ...
                        strcmp(antenna(istat).name, station_name_old(ista)))
                            delta_ZWD(icolo) = delta_ZWD_f(ista);  % 1*ZWD_Hb - 1*ZWD_Ho = Delta D
                    end
                end
            end
        end
    end
    if antenna(istat).name == 'HARTRAO '
        for jstat = 1:na
            if antenna(jstat).name == 'HART15M '
                icolo=icolo+1;
                ncombi(icolo,1)=istat;
                ncombi(icolo,2)=jstat;
                for ista=1:nsta 
                    if (strcmp(antenna(jstat).name, station_name_new(ista)) && ...
                        strcmp(antenna(istat).name, station_name_old(ista)))
                            delta_ZWD(icolo) = delta_ZWD_f(ista);
                    end
                end
            end
        end
    end
    if antenna(istat).name == 'WETTZELL'
        for jstat = 1:na
            if antenna(jstat).name == 'WETTZ13N'
                icolo=icolo+1;
                ncombi(icolo,1)=istat;
                ncombi(icolo,2)=jstat;
                for ista=1:nsta 
                    if (strcmp(antenna(jstat).name, station_name_new(ista)) && ...
                        strcmp(antenna(istat).name, station_name_old(ista)))
                            delta_ZWD(icolo) = delta_ZWD_f(ista);
                    end
                end
            end
        end
    end
    if antenna(istat).name == 'YEBES40M'
        for jstat = 1:na
            if antenna(jstat).name == 'RAEGYEB '
                icolo=icolo+1;
                ncombi(icolo,1)=istat;
                ncombi(icolo,2)=jstat;
                for ista=1:nsta 
                    if (strcmp(antenna(jstat).name, station_name_new(ista)) && ...
                        strcmp(antenna(istat).name, station_name_old(ista)))
                            delta_ZWD(icolo) = delta_ZWD_f(ista);
                    end
                end
            end
        end
    end
end 

%-----------------------------------------%

Hcomzwd = []; Phcomzwd = []; ochcomzwd = [];
t_.zwd = sum([n_.zwd]); % total number of zwd paramters in a session
ncolo=icolo; % number of co-location sites % currently only 4, Hobart, Hartebeesthoek, Wettzell, Yebes
sigma=33.3e-12; % sigma in weight matrix in seconds e.g. 33ps (around 1cm)  
for icom = 1:ncolo % number of co-location site
    sumzwdj = 0;
    sumzwdk = 0;

    j=ncombi(icom,1);  % istat for co-location site  !! the order at antenna.mat HOBART26(Ho)
    k=ncombi(icom,2);  % istat for co-location site  !! the order at antenna.mat HOBART12(Hb)
    H_comzwd(icom).h(n_(j).zwd,t_.zwd) = 0;     % initialization for constraint matrix        !! the matrix size: (1) row    = no. of co-location sites * number of zwd per
                                                %                                                                 (2) column = no. of zwd parameters 
    %Ph_comzwd(icom).h(n_(j).zwd,n_(j).zwd)= 0;  % initialization for constraint weight matrix !! diagonal square matrix with same number of row of H_zwd               
    %oc_hcomzwd(icom).h(n_(j).zwd,1) = 0;

    for i = 1:j-1 % for both antennas j & k
        sumzwdj = sumzwdj + n_(i).zwd;  % counting number of zwd estimates right before the estimate to use
                                        % if it is the 6th site, count zwd numbers till end of 5th site
    end
    for i = 1:k-1 
        sumzwdk = sumzwdk + n_(i).zwd;
    end

    inter=0;
    for jnter = 1:n_(j).zwd
        for knter = 1:n_(k).zwd
            if (t(j).zwd(jnter) == t(k).zwd(knter))
                inter=inter+1;
                H_comzwd(icom).h(inter,sumzwdj+jnter) = -1;  %  e.g. zwd of 'HOBART26(Ho)' % 1*ZWD_Hb - 1*ZWD_Ho = Delta D    
                H_comzwd(icom).h(inter,sumzwdk+knter) = +1;  %  e.g. zwd of 'HOBART12(Hb)' 
                Ph_comzwd(icom).h(inter,inter) = 1./...
                    ((sigma)^2*c^2*100^2); % [1/cm^2] weight matrix of the design matrix for the zwd tie constarints
                oc_hcomzwd(icom).h(inter,1) = delta_ZWD(icom);  % ZWD difference according to height  
            end
        end
    end
    % Concatenating
    Hcomzwd = vertcat(Hcomzwd,H_comzwd(icom).h);
    Phcomzwd = blkdiag(Phcomzwd,Ph_comzwd(icom).h);
    ochcomzwd = vertcat(ochcomzwd,oc_hcomzwd(icom).h);
end

%save h.mat H_comzwd   % test
%save ph.mat Ph_comzwd %test

% add to original constrain matrix
H(3).sm   = vertcat(H(3).sm,   Hcomzwd);   
Ph(3).sm  = blkdiag(Ph(3).sm,  Phcomzwd); 
och(3).sv = vertcat(och(3).sv, ochcomzwd); 
