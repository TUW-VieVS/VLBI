% ************************************************************************
% Description:
% This function corrects the effect of galactic aberration to the source coordinates.
% A similar approach is used by nuSolve
%
%
%   Input:
%      ---
%      the program is called from vie_mod.m
%      Data input:
%       - list of RA 
%       - list of DEC  
%       - mean MJD of session
%
%   Output:
%       - list of RA 
%       - list of DEC  
%
%
%   Coded for VieVS:
%   06 Sep 2019 by David Mayer


function [de_GA, ra_GA] = correct_GA(de, ra, ave_date_session)

    RA_Gal = deg2rad(266.4); % coordinates of galactic center
    DE_Gal = deg2rad(-28.94);
    GA_val = 5.8; %muas/year
    GA_ref = date2mjd([2015 1 1]); %reference epoch of icrf3
    GA_val = deg2rad(GA_val*1e-6/3600); %rad/year
    
    GA_vec = [  cos(DE_Gal)*cos(RA_Gal);
                cos(DE_Gal)*sin(RA_Gal);
                sin(DE_Gal)]; % coordinates of galactic center
    time_delta = (ave_date_session - GA_ref)/365.25; % differrence to reference epoch in years

    SIN_DE = sin(de);
    COS_DE = cos(de);
    SIN_RA = sin(ra);
    COS_RA = cos(ra);  
    
    delta_RA = (-GA_vec(1).*SIN_RA          + GA_vec(2).*COS_RA)./COS_DE;
    delta_DE =  -GA_vec(1).*SIN_DE.*COS_RA  - GA_vec(2).*SIN_DE.*SIN_RA     + GA_vec(3).*COS_DE;
    
    delta_RA = delta_RA.*GA_val; 
    delta_DE = delta_DE.*GA_val;
   
    ra_GA = ra + delta_RA*time_delta;
    de_GA = de + delta_DE*time_delta;
