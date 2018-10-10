%**************************************************************************
%
% Program to calculate the variance per baseline
%
% Input: 
%   - scan
%   - v_real
%   - antenna
%
% Output:
%   - var_bas   (1,number of baselines)
%   - Pobserv
%
% Coded for VieVS (adapted from calc_chi2.m):
% Minttu Uunila (8 Jul 2014)
%
% Changes:
% 9 Apr 2015 by Hana Krasna: bugs in indexing corrected
%
%**************************************************************************

function [var_bas, Pobserv] = calc_bas_weights(scan,v_real,antenna)



% fid=fopen(['../TEMP/Bslweights/VieVS_baseline_weights/' antenna(1).session(1:9) '_v004.txt'],'wt');
% fprintf(fid,'%s   \n' ,'! Baseline              G-Delay   G-Rate   P-Delay    P-Rate ');
% 


c = 299792458;
% number of scans
ntim = length(scan); 
% number of antennas
na = length(antenna);
nobserv = sum([scan.nobs]);
temp = [scan.obs];
mi_observ = [temp.sig]; % [seconds]
temp =[];


% Creates matrix stat containing the numbers of the stations participating
% in every observation

stat(1:ntim,1:2)=0; % Preallocating
i=1;                % Initializing counter
for nscan=1:ntim
   for nob=1:scan(nscan).nobs
      stat(i,1)=scan(nscan).obs(nob).i1;
      stat(i,2)=scan(nscan).obs(nob).i2;
      i=i+1; 
   end
end
            
disp(' ')

% Initializing of the P matrix with baseline INDEPENDENT values
addnoise=0.01;
nmi_observ = sqrt((addnoise/c).^2.+(mi_observ).^2); % [seconds]
Pobserv = diag(sparse(1./((nmi_observ.^2).*c^2*100^2))); % [1/cm^2]


for i1=1:na-1
    for i2=i1+1:na
        v=[];
        ids=[]; ki=0;
        for j=1:nobserv
            if (stat(j,1)==i1 || stat(j,2)==i1) && (stat(j,1)==i2 || stat(j,2)==i2)
                v = [v;v_real(j,1)];
                ki=ki+1;
                ids(ki)=j;
            end
        end
   
        if ~isempty(v)
            n=length(v);
            var_bas=[];
            var_bas=(sum(v.^2)/n);

           disp(sprintf('variance of baseline %s - %s cm: %3.4f ',antenna(i1).name,antenna(i2).name, sqrt(var_bas)));
%              disp(sprintf('variance of baseline %s - %s ps: %3.4f ',antenna(i1).name,antenna(i2).name, sqrt(var_bas)*(0.01/c)*10^12));
            
%%
%             idblank=setdiff(1:8,1:max(find(antenna(i1).name(1,:)~=' ')));
%             a1=antenna(i1).name;
%             a1(idblank)='_';
%            
%             idblank=setdiff(1:8,1:max(find(antenna(i2).name(1,:)~=' ')));
%             a2=antenna(i2).name;
%             a2(idblank)='_';
%             
% 
%             fprintf(fid,'  %c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c     %6.2f    %6.2f    %6.2f    %6.2f   \n' ,a1,'/',a2,sqrt(var_bas)*(0.01/c)*10^12,0,0,0);
            
            
%%

          

            for k=1:ki
                 nmi_observ=sqrt(var_bas*(0.01/c)^2+mi_observ(ids(k)).^2); %s
                 Pobserv(ids(k),ids(k))=1./((nmi_observ.^2).*c^2*100^2);
            end
            clear n v 
        end
    end
end

% fclose(fid);