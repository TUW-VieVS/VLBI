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
%
%   Output:
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   12 May 2009 by Kamil Teke
%
%   Revision: 
%   06 Dec 2009 by Kamil Teke: header added
%   10 Aug 2012 by Kamil Teke: units of constraints are now in mas & cm
%   10 Oct 2016 by A. Girdiuk: code-optimized
% ************************************************************************
function [H,Ph,och] = hpoc_xyz(H,Ph,och,n_,opt,na)

Hx = []; Phx = []; Hy = []; Phy = []; Hz = []; Phz = [];
%t_.xyz = sum([n_.xyz]); % total estimate interval of coordinates in a session


% FORMING THE CONSTRAINTS of COORDINATES
if  opt.pw_stc == 0
    Hx(na-1,na) = 0; %Hy = Hx; Hz = Hx;
    Phx(na-1,na-1) = 0;% Phy = Phx; Phz = Phx;   
elseif opt.pw_stc == 1
    for istat = 1:na  
        coef_xyz = opt.stat(istat).coef_xyz;
        
        mat_xyz = diag(ones(1,n_(istat).xyz)) - diag(ones(1,n_(istat).xyz-1),1);
        Hx= blkdiag(Hx,mat_xyz(1:n_(istat).xyz-1,1:n_(istat).xyz));
        Phx = blkdiag(Phx,diag(ones(1,n_(istat).xyz-1).*1./coef_xyz^2));
    end
end

Hy = Hx; 
Phy = Phx; 
Hz = Hx; 
Phz = Phx;

% -------------------------------------------------------------------------
% FORMING THE O-C VECTOR FOR THE CONSTRAINTS
oc_hx(size(Hx,1),1) = 0;  % O-C vector for the constraints of  x coordinates
oc_hy(size(Hy,1),1) = 0;  % O-C vector for the constraints of  y coordinates
oc_hz(size(Hz,1),1) = 0;  % O-C vector for the constraints of  z coordinates

if opt.constr_xyz == 0
    Hx(1:size(Hx,1),1:size(Hx,2)) = 0; Hy = Hx; Hz = Hx; 
    Phx(1:size(Phx,1),1:size(Phx,2)) = 0; Phy = Phx; Phz = Phx;
    oc_hx = []; oc_hy = []; oc_hz = [];
end

H(13).sm = Hx; H(14).sm = Hy; H(15).sm = Hz;
Ph(13).sm = Phx; Ph(14).sm = Phy; Ph(15).sm = Phz;

och(13).sv = oc_hx; och(14).sv = oc_hy; och(15).sv = oc_hz;
