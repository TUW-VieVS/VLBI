% function to calculate the contribution of source structure effects
% Lucia, 30-04-2014
% Lucia, 12-09-2014: IERS/IVS source name

function [obs_ok,grout]= calc_soustruc_pscan(scan,source,station)

clear jetangpsc
% obs_ok = 0 ...... this scan is bad
% obs_ok = 1 ...... this scan is good


% SOURCE STRUCTURE +
 
% read in source structure catalogue 
    % filename_sou_cat=strcat('../CRF/SOURCE_STRUCTURE_CAT/sim_struct_2comp_fluxRatio_random_ICRF2');
%     filename_sou_cat=strcat('D:\VieVS_jetsched\CRF\SOURCE_STRUCTURE_CAT\sim_struct_2comp_fluxRatio_random_ICRF2.cat');
 filename_sou_cat=strcat('../CRF/SOURCE_STRUCTURE_CAT/sim_struct_2comp_fluxRatio_random_SI_3.cat');
    [cat_comp.name,cat_comp.flux,cat_comp.maj,cat_comp.min,cat_comp.angle,cat_comp.dRA,cat_comp.dDec]=textread(filename_sou_cat,'%s%f%f%f%f%f%f','delimiter',',');    
 
% loop over all observations in a scan
nst = scan.nsta;   %%%
nobs=(nst*(nst-1))/2;
obs=1;
for io=1:nst-1   %%%
    for iio=(io+1):nst   %%%
        st1 = scan.sta(io).staid; % index of 1st station   %%%
        st2 = scan.sta(iio).staid;   %%%
        ind = strmatch(source(scan.srcid).name,cat_comp.name,'exact');
            if isempty(ind)
                nam2=sounamivs2iers(source(scan.srcid).name);
                ind = strmatch(nam2,cat_comp.name,'exact');
            end
            if isempty(ind)
                disp(strcat('source ',source(scan.srcid).name,' not found in ss catalogue.'));   %%%
                soucorr=0; jetang=90; jetjb=0; uvrange=0; uu=0; vv=0;
            else 
                source(scan.srcid).sou_model=[cat_comp.flux(ind),cat_comp.maj(ind),cat_comp.min(ind),cat_comp.angle(ind),cat_comp.dRA(ind),cat_comp.dDec(ind)];
                [soucorr,uu,vv]=modDelay(source(scan.srcid).sou_model,station(st1).xyz,station(st2).xyz,([8213 8252 8353 8513 8733 8853 8913 8933]+4), 8217, source(scan.srcid).ra, source(scan.srcid).de, scan.startmjd);   %%%
                jetvec=(source(scan.srcid).sou_model(2,5:6));   %%%
                jetvec=jetvec/norm(jetvec);
                uvvec=[uu;vv];
                uvvec=uvvec/norm(uvvec);          
                jetang=acos(abs(jetvec*uvvec))*180/pi;  
                uvrange=(jetvec)*[uu;vv];                
            end
    jetangpsc(obs)=jetang;
    uvrangpsc(obs)=uvrange;
    obs=obs+1;
    end
end

 gr=mean(jetangpsc);
 grout=gr/90;
%gr=median(abs(uvrangpsc));

%grout=gr/350;

% gr=sum(jetangpsc)/nobs;
%   meduv(ksc)=sum(abs([scan(isc).obs.uvrange]))/(scan(isc).nobs);


if gr > 40
%if gr < 100    
    obs_ok=1;
else
    obs_ok=0;
end

