% ************************************************************************
%   Description:
%	Reads the scan struct and organises it to observations per station
%
%   Input:										
%     scan              scan struct
%     na                number of antennas
%     opt				settings
%     mjd0				beginning of the session
%     antenna			antenna struct
%
% 
%   Output:
%     per_stat          struct including data per station wise
%     obs_mjd           time of observation
% 
%   External calls: 	
%       
%   Coded for VieVS: 
%   17 Dec 2024 by Helene Wolf
%
%   Revision: 
%
%
% ************************************************************************
function [per_stat, obs_mjd] = scan2obsPerStation(scan, na, opt, mjd0, antenna)
    % preallocation
    if ~opt.est_stsespos
        per_stat(1,na)=struct('mjd', [], 'oc_nob', [], 'zd', [], 'first', [], 'other', [], 'mf', [], 'az', [],...
                     'dx', [],     'dy', [], 'dz', [],    'xo', [],    'yo', [], 'zo', [],...
                     'dAO', [], 'drg', []);
    else
        per_stat(1,na)=struct('mjd', [], 'oc_nob', [], 'zd', [], 'first', [], 'other', [], 'mf', [], 'az', [],...
                     'dx', [],     'dy', [], 'dz', [],    'xo', [],    'yo', [], 'zo', [],...
                     'dAO', [], 'drg', [], 'pAcr', [], 'pAce', [], 'pAcn', [], 'pAsr', [], 'pAse', [], 'pAsn', []);
    end
    
    n_scan = length(scan);
    for istat = 1:na 
        i = 0; n_obs_per_src = 0; n_obs_per_sat=0; n_obs_per_qu = 0;
        for itim = 1:n_scan 
            for iobs = 1:scan(itim).nobs % number of observations per scan
                i = i + 1; % i : i. observation in the session
                i1 = scan(itim).obs(iobs).i1;
                i2 = scan(itim).obs(iobs).i2;
                if i1 == istat || i2 == istat
                    n_obs_per_src = n_obs_per_src + 1; % k : k. observation of the specific station
                    per_stat(istat).mjd(n_obs_per_src) = scan(itim).mjd; % The times of scans [day]
                    per_stat(istat).oc_nob(n_obs_per_src) = i;
                    per_stat(istat).zd(n_obs_per_src) = scan(itim).stat(istat).zd;  % Boehm 21 Aug 2009, 15 Sep 2010
                    obs_mjd(i) = (scan(itim).mjd - mjd0)*24*60; % [minute]
                    per_stat(istat).obs_type(n_obs_per_src)  = string(scan(itim).obs_type);
                    if strcmp (scan(itim).obs_type,"s")
                       n_obs_per_sat = n_obs_per_sat + 1;
                       per_stat(istat).oc_nob_sat(n_obs_per_sat) = i;
                    else
                       n_obs_per_qu = n_obs_per_qu + 1;
                       per_stat(istat).oc_nob_qu(n_obs_per_qu) = i;
                    end
                end
                if i1 == istat
                    per_stat(istat).first(n_obs_per_src) = -1;
                    per_stat(istat).other(n_obs_per_src) = i2;  % Boehm 21 Aug 2009
    
                    if strcmp(scan(itim).obs_type, "s")
                        per_stat(istat).first_sat(n_obs_per_sat)  = -1;
                        per_stat(istat).dx_sat(n_obs_per_sat) = -scan(itim).obs(iobs).pstat1(1); % The partial derivatives of delay wrt antenna coordinates
                        per_stat(istat).dy_sat(n_obs_per_sat) = -scan(itim).obs(iobs).pstat1(2);
                        per_stat(istat).dz_sat(n_obs_per_sat) = -scan(itim).obs(iobs).pstat1(3);
                    else
                        per_stat(istat).first_qu(n_obs_per_qu)  = -1;
                        per_stat(istat).dx_qu(n_obs_per_qu) = -scan(itim).obs(iobs).pstat1(1); % The partial derivatives of delay wrt antenna coordinates
                        per_stat(istat).dy_qu(n_obs_per_qu) = -scan(itim).obs(iobs).pstat1(2);
                        per_stat(istat).dz_qu(n_obs_per_qu) = -scan(itim).obs(iobs).pstat1(3);
                    end
    
                    per_stat(istat).mf(n_obs_per_src) = scan(itim).stat(i1).mfw; % The mapping function value
                    per_stat(istat).az(n_obs_per_src) = scan(itim).stat(i1).az; % Azimuth [radians]
                    %per_stat(istat).zd(k) = scan(itim).stat(i1).zd; % Zenith distance [radians]
    
                    per_stat(istat).dx(n_obs_per_src) = -scan(itim).obs(iobs).pstat1(1); % The partial derivatives of delay wrt antenna coordinates
                    per_stat(istat).dy(n_obs_per_src) = -scan(itim).obs(iobs).pstat1(2);
                    per_stat(istat).dz(n_obs_per_src) = -scan(itim).obs(iobs).pstat1(3);
    
                    %obs_per_stat.dx(k) = scan(itim).stat(i1).pantd(1); % The partial derivatives of delay wrt antenna coordinates
                    %obs_per_stat.dy(k) = scan(itim).stat(i1).pantd(2);
                    %obs_per_stat.dz(k) = scan(itim).stat(i1).pantd(3);
    
                    per_stat(istat).xo(n_obs_per_src) = scan(itim).stat(i1).x(1); % Apriori coordinates of the antennas
                    per_stat(istat).yo(n_obs_per_src) = scan(itim).stat(i1).x(2); % y
                    per_stat(istat).zo(n_obs_per_src) = scan(itim).stat(i1).x(3); % z
    
                    per_stat(istat).dAO(n_obs_per_src) = -scan(itim).obs(iobs).pAO_st1;
    
                    per_stat(istat).drg(n_obs_per_src)= -scan(itim).obs(iobs).prg_st1;
    
                     if opt.est_stsespos ==1
    		            per_stat(istat).pAcr(n_obs_per_src,:) = -scan(itim).obs(iobs).pAcr_st1;
    		            per_stat(istat).pAce(n_obs_per_src,:) = -scan(itim).obs(iobs).pAce_st1;
    		            per_stat(istat).pAcn(n_obs_per_src,:) = -scan(itim).obs(iobs).pAcn_st1;
    		            per_stat(istat).pAsr(n_obs_per_src,:) = -scan(itim).obs(iobs).pAsr_st1;
    		            per_stat(istat).pAse(n_obs_per_src,:) = -scan(itim).obs(iobs).pAse_st1;
		                per_stat(istat).pAsn(n_obs_per_src,:) = -scan(itim).obs(iobs).pAsn_st1;
                     end
                end
                if i2 == istat
                    per_stat(istat).first(n_obs_per_src) = +1;
    
                    if strcmp(scan(itim).obs_type, "s")
                        per_stat(istat).first_sat(n_obs_per_sat)  = +1;
                        per_stat(istat).dx_sat(n_obs_per_sat) = scan(itim).obs(iobs).pstat2(1); % The partial derivatives of delay wrt antenna coordinates
                        per_stat(istat).dy_sat(n_obs_per_sat) = scan(itim).obs(iobs).pstat2(2);
                        per_stat(istat).dz_sat(n_obs_per_sat) = scan(itim).obs(iobs).pstat2(3);
                    else
                        per_stat(istat).first_qu(n_obs_per_qu)  = +1;
                        per_stat(istat).dx_qu(n_obs_per_qu) = scan(itim).obs(iobs).pstat2(1); % The partial derivatives of delay wrt antenna coordinates
                        per_stat(istat).dy_qu(n_obs_per_qu) = scan(itim).obs(iobs).pstat2(2);
                        per_stat(istat).dz_qu(n_obs_per_qu) = scan(itim).obs(iobs).pstat2(3);   
                    end
                    per_stat(istat).other(n_obs_per_src) = i1;  % Boehm 21 Aug 2009
    
                    per_stat(istat).mf(n_obs_per_src) = scan(itim).stat(i2).mfw; % The mapping function value
                    per_stat(istat).az(n_obs_per_src) = scan(itim).stat(i2).az; % Azimuth [radians]
    %                    per_stat(istat).zd(k) = scan(itim).stat(i2).zd; % Zenith distance [radians]
    
                    per_stat(istat).dx(n_obs_per_src) = scan(itim).obs(iobs).pstat2(1); % The partial derivatives of delay wrt antenna coordinates
                    per_stat(istat).dy(n_obs_per_src) = scan(itim).obs(iobs).pstat2(2);
                    per_stat(istat).dz(n_obs_per_src) = scan(itim).obs(iobs).pstat2(3);
    
                    %obs_per_stat.dx(k) = scan(itim).stat(i2).pantd(1); % The partial derivatives of delay wrt antenna coordinates
                    %obs_per_stat.dy(k) = scan(itim).stat(i2).pantd(2);
                    %obs_per_stat.dz(k) = scan(itim).stat(i2).pantd(3);
    
                    per_stat(istat).xo(n_obs_per_src) = scan(itim).stat(i2).x(1); % Apriori coordinates of the antennas
                    per_stat(istat).yo(n_obs_per_src) = scan(itim).stat(i2).x(2);
                    per_stat(istat).zo(n_obs_per_src) = scan(itim).stat(i2).x(3);
    
                    per_stat(istat).dAO(n_obs_per_src) = scan(itim).obs(iobs).pAO_st2;
    
                    per_stat(istat).drg(n_obs_per_src) = scan(itim).obs(iobs).prg_st2;
    
                     if opt.est_stsespos ==1
        		        per_stat(istat).pAcr(n_obs_per_src,:) = scan(itim).obs(iobs).pAcr_st2;
        		        per_stat(istat).pAce(n_obs_per_src,:) = scan(itim).obs(iobs).pAce_st2;
        		        per_stat(istat).pAcn(n_obs_per_src,:) = scan(itim).obs(iobs).pAcn_st2;
        		        per_stat(istat).pAsr(n_obs_per_src,:) = scan(itim).obs(iobs).pAsr_st2;
        		        per_stat(istat).pAse(n_obs_per_src,:) = scan(itim).obs(iobs).pAse_st2;
		                per_stat(istat).pAsn(n_obs_per_src,:) = scan(itim).obs(iobs).pAsn_st2;
                     end
                end
            end
        end
        fprintf('obs. of the antenna %s : %4d\n',antenna(istat).name,n_obs_per_src);
    end
end