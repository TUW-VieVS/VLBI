% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   constraints of antennas' coordinates
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
%   12 May 2009 by Kamil Teke
%
%   Revision: 
%   06 Dec 2009 by Kamil Teke: header added
%   10 Aug 2012 by Kamil Teke: units of constraints are now in mas & cm
%   28 May 2015 by Lucia Plank: hardcoded option of constraining HbHo
% ************************************************************************
function [H,Ph,och,t_] = hpoc_xyz_hbho(H,Ph,och,n_,opt,na,t_,HbHo)

Hx = []; Phx = []; Hy = []; Phy = []; Hz = []; Phz = [];
t_.xyz = sum([n_.xyz]); % total estimate interval of coordinates in a session

fprintf('INTRODUCE CONSTRAINTS FOR HOBART12-HOBART26\n');
gewicht = 1; % [cm] % weight


for istat = 1:na  
    % FORMING THE CONSTRAINTS of COORDINATES
    if  opt.pw_stc == 0
        Hx(na-1+1,na) = 0; 
        Hx(na-1+1,HbHo(1))=1;   Hx(na-1+1,HbHo(2))=-1;
        Hy = Hx; Hz = Hx;
%         Phx(na-1+1,na-1) = 0; 
        Phx(na,na) = 0; 
        Phx(na,na)=1/gewicht^2;
%         Phx(na,HbHo(1))=1/gewicht^2; Phx(na,HbHo(2))=1/gewicht^2;
        Phy = Phx; Phz = Phx;   
    else if opt.pw_stc == 1
    coef_xyz = opt.stat(istat).coef_xyz;
    int_xyz = opt.stat(istat).int_xyz;
    H_x(istat).h(t_.xyz-na,n_(istat).xyz) = 0; 
    P_hx(istat).h(t_.xyz-na,n_(istat).xyz-1)= 0;        
    sumxyz = 0;
        for i = 1:istat-1        
            sumxyz = sumxyz + n_(i).xyz-1;
        end      
        for inter = 1:n_(istat).xyz-1 
            H_x(istat).h(sumxyz + inter, inter) = +1;  % design matrix for the x coordinates pseudo observtaions as constraints
            H_x(istat).h(sumxyz + inter, inter+1) = -1; % design matrix for the x coordinates pseudo observtaions as constraints
            P_hx(istat).h(sumxyz + inter, inter) = 1./coef_xyz^2; % [1/cm^2] weight matrix of the design matrix for the coordinates constraints
        end          
        Hx = horzcat(Hx,H_x(istat).h); % Concatenating
        Phx = horzcat(Phx,P_hx(istat).h); % Concatenating
        Hy = Hx; 
        Phy = Phx; 
        Hz = Hx; 
        Phz = Phx;
        end
    end
end

% -------------------------------------------------------------------------
% FORMING THE O-C VECTOR FOR THE CONSTRAINTS
oc_hx(size(Hx,1),1) = 0;  % O-C vector for the constraints of  x coordinates
oc_hy(size(Hy,1),1) = 0;  % O-C vector for the constraints of  y coordinates
oc_hz(size(Hz,1),1) = 0;  % O-C vector for the constraints of  z coordinates

% oc_hx(end)=2465.743;
% oc_hy(end)= 735.364;
% oc_hz(end)=1461.538;

oc_hx(end)=0;
oc_hy(end)= 0;
oc_hz(end)=0;

if opt.constr_xyz == 0
    Hx(1:size(Hx,1),1:size(Hx,2)) = 0; Hy = Hx; Hz = Hx; 
    Phx(1:size(Phx,1),1:size(Phx,2)) = 0; Phy = Phx; Phz = Phx;
    oc_hx = []; oc_hy = []; oc_hz = [];
end

H(13).sm = Hx; H(14).sm = Hy; H(15).sm = Hz;
Ph(13).sm = Phx; Ph(14).sm = Phy; Ph(15).sm = Phz;

och(13).sv = oc_hx; och(14).sv = oc_hy; och(15).sv = oc_hz;