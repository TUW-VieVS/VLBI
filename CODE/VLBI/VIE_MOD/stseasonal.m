% ************************************************************************
%   Description:
%   Amplitudes of seasonal variations in
%   station positions (annual and semi-annual)
%  
%   Input:										
%      mjd                 Modified Julian Date [d]
%      phi,lam             lat, long of the stations
% 
%   Output:
%      partials w.r.t. the cosine and sine amplitudes of the station motion in xyz
% 
%
%   Coded for VieVS: 
%   11 Oct 2011 by Hana Spicakova: 


function [pAcr_xyz pAce_xyz pAcn_xyz pAsr_xyz pAse_xyz pAsn_xyz]= stseasonal(mjd,phi,lam)

mjd0=51544; %2000

%cpsd
freq_cpsd= [
     0.005461        %Ssa - semi-annual
     0.002730375267416];   %365.25 solar days - annual


nt = length(freq_cpsd);
 
sol2sid=366.2422/365.2422; %1.002737909

Psd = 1./freq_cpsd; % sid day
P = Psd./sol2sid; % solar day

% zero apriori
% dr = Ac_r.*cos((mjd-mjd0)./P*2*pi) + As_r.*sin((mjd-mjd0)./P*2*pi);
% de = Ac_e.*cos((mjd-mjd0)./P*2*pi) + As_e.*sin((mjd-mjd0)./P*2*pi);
% dn = Ac_n.*cos((mjd-mjd0)./P*2*pi) + As_n.*sin((mjd-mjd0)./P*2*pi);

% partials for cosinus amplitudes
pAcr = [cos((mjd-mjd0)./P*2*pi) zeros(nt,1) zeros(nt,1) ];
pAce = [zeros(nt,1) cos((mjd-mjd0)./P*2*pi) zeros(nt,1)];
pAcn = [zeros(nt,1) zeros(nt,1)             cos((mjd-mjd0)./P*2*pi)];

% partials for sinus amplitudes
pAsr = [sin((mjd-mjd0)./P*2*pi) zeros(nt,1) zeros(nt,1) ];
pAse = [zeros(nt,1) sin((mjd-mjd0)./P*2*pi) zeros(nt,1)];
pAsn = [zeros(nt,1) zeros(nt,1)             sin((mjd-mjd0)./P*2*pi)];


% transform each amplitude separately into xyz
phint(1:nt,1)=phi;
lamnt(1:nt,1)=lam;
pAcr_xyz =ren2xyz(pAcr,phint,lamnt);
pAce_xyz =ren2xyz(pAce,phint,lamnt);
pAcn_xyz =ren2xyz(pAcn,phint,lamnt);

pAsr_xyz =ren2xyz(pAsr,phint,lamnt);
pAse_xyz =ren2xyz(pAse,phint,lamnt);
pAsn_xyz =ren2xyz(pAsn,phint,lamnt);

