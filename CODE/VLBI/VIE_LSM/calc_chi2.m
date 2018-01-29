%**************************************************************************
%
% Program to calculate the chi2 per station and per baseline
%
% Input: 
%   - na
%   - ntim
%   - nobserv
%   - scan
%   - Pobserv
%   - v_real
%   - antenna
%
% Output:
%   - mo_stat  (1,na) 
%   - mo_bas   (1,number of baselines)
%
% Coded for VieVS:
% Claudia Tierno Ros (13 Nov. 2012)
%
%   Revision: 
%   26 Apr 2013 by Sigrid Boehm: screen output of m0 changed to chi^2
%**************************************************************************

function [mo_stat mo_bas]=calc_chi2(na,ntim,nobserv,scan,Pobserv,v_real,antenna)

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

% chi2 per station
mo_stat(1,na)=0;
for i=1:na
 v=[];
 P=[];  

    for j=1:nobserv
        if stat(j,1)==i || stat(j,2)==i
            v = [v;v_real(j,1)];
            P = blkdiag(P,Pobserv(j,j));
        end
    end
 
    n=length(v);
    vTPv=v'*P*v;
    mo_stat(1,i)=sqrt(vTPv/n);
    clear n vTPv v P
    
    disp(sprintf('chi-squared of antenna %s : %3.4f ',antenna(i).name, mo_stat(1,i).^2));
    
end
            
disp(' ')

% chi2 per baseline
nbas=(na*(na-1))/2;
mo_bas(1,1:nbas)=999;
ibas=1;

for i1=1:na-1
    
    for i2=i1+1:na
        v=[];
        P=[];
        for j=1:nobserv
            if (stat(j,1)==i1 || stat(j,2)==i1) && (stat(j,1)==i2 || stat(j,2)==i2)
                v = [v;v_real(j,1)];
                P = blkdiag(P,Pobserv(j,j));
            end
        end
        
        if ~isempty(v)
        n=length(v);
        vTPv=v'*P*v;
        mo_bas(1,ibas)=sqrt(vTPv/n);
        clear n vTPv v P
        disp(sprintf('chi-squared of baseline %s - %s: %3.4f ',antenna(i1).name,antenna(i2).name, mo_bas(1,ibas).^2));
        end
          
        ibas=ibas+1;  
    end
end


