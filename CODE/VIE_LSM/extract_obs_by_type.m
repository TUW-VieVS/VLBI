% ************************************************************************
%   Description:
%	get the station coordinates in East, North and Up directions
%
%
%   Input:										
%     scan                scan structure
%     n_sources           number of sources for this type
%     type_char           type of source (s for satellite and q for quasar)
% 
% 
%   Output:
%     obs_per_item          observations sorted per item
% 
%   External calls: 	
%
%       
%   Coded for VieVS (taken from vie_lsm): 
%   10 June 2026 by Helene Wolf
%
%   Revision: 
%
% ************************************************************************

function obs_per_item = extract_obs_by_type(scan, n_sources, type_char)
    % get relevant scans
    valid_mask = strcmp({scan.obs_type}, type_char);
    rel_scans = scan(valid_mask);
    
    if isempty(rel_scans)
        obs_per_item = [];
        return;
    end
    
    unique_ids = unique([rel_scans.iso]);
    if ~isempty(n_sources) && isnumeric(n_sources)
        unique_ids = intersect(unique_ids, 1:n_sources);
    end
    
    n_src = numel(unique_ids);
    obs_per_item = repmat(struct('mjd',[],'iso',[],'i1',[],'i2',[]), 1, n_src);
    
    for k = 1:n_src
        src_id = unique_ids(k);
        
        idx = [rel_scans.iso]' == repmat(src_id, length(rel_scans),1);
        src_scans = rel_scans(idx);
        
        n_total = sum([src_scans.nobs]);
        if n_total == 0, continue; end
        
        obs_per_item(k).mjd = zeros(1, n_total);
        obs_per_item(k).iso = zeros(1, n_total);
        obs_per_item(k).i1  = zeros(1, n_total);
        obs_per_item(k).i2  = zeros(1, n_total);
        
        pos = 1;
        for j = 1:numel(src_scans)
            len = src_scans(j).nobs;
            obs_per_item(k).mjd(pos:pos+len-1) = src_scans(j).mjd;
            obs_per_item(k).iso(pos:pos+len-1) = src_scans(j).iso;
            obs_per_item(k).i1(pos:pos+len-1)  = [src_scans(j).obs(1:len).i1];
            obs_per_item(k).i2(pos:pos+len-1)  = [src_scans(j).obs(1:len).i2];
            pos = pos + len;
        end
    end
end