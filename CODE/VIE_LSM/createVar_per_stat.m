function [per_stat, obs_mjd] = createVar_per_stat(scan, antenna, opt, per_stat, na, n_scan, mjd0, type)

for istat = 1:na % number of stations
    i = 0; n_obs_per_src = 0;
    for itim = 1:n_scan % number of scans per session       
        %if strcmp(scan(itim).obs_type, type) || strcmp(type,'all')
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
                    per_stat(istat).obs_type = scan(itim).obs_type;
                end
                if i1 == istat
                    per_stat(istat).first(n_obs_per_src) = -1;
                    per_stat(istat).other(n_obs_per_src) = i2;  % Boehm 21 Aug 2009
    
                    %
                    %   partials prepared
                    %
    
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
           % end
        end
    end
end
 if ~strcmp(scan(itim).obs_type, type)
    fprintf('obs. of the antenna %s : %4d\n',antenna(istat).name,n_obs_per_src);
 end

end