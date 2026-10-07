% ************************************************************************
%   Description:
%   Turbulence-based Variance-Covariance Matrix for VLBI observations 
%
%   Input:
%      parameter  - VieVS parameter struct
%      mi_observ  - formal observation accuracy [s]
%      cLight     - speed of light [m/s]
%      antenna    - antenna struct 
%      sources    - source struct 
%      scan       - scan struct with observations, elevations, azimuths
%
%   Output:
%      addnoise       - placeholder scalar (= 1)
%      Sigma_turb_all - full Variance-Covariance Matrix [cm^2]
%
%   Reference:
%      Halsig, S. (2018). Atmospheric refraction and turbulence in VLBI data analysis (Doctoral 
%      dissertation, Rheinische Friedrich-Wilhelms-Universität Bonn, Landwirtschaftliche Fakultät, 
%      IGG-Institut für Geodäsie und Geoinformation).
%      Kermarrec, G., & Schön, S. (2014). On the Mátern covariance family: a proposal for modeling 
%      temporal correlations based on turbulence theory. Journal of Geodesy, 88(11), 1061-1079.
%
%   Coded for VieVS:
%       07.10.2026 by Peter Urban
%
%   Update:
%       
%
% ************************************************************************

function [addnoise, Sigma_turb_all] = vievs_turb_model(parameter, mi_observ, cLight, antenna, sources, scan)

tic;
format longg

% Read and store all observations 
totalObservations = sum(arrayfun(@(x) length(x.obs), scan));
all_obs(totalObservations) = struct('orig_idx', [], 'scan', [], 'obs', [], 'tau', [], 'tau_sig', [], ...
    'i1', [], 'i2', [], 'mjd', [], 'e1', [], 'e2', [], 'az1', [], 'az2', [], ...
    'x1', [], 'y1', [], 'z1', [], 'x2', [], 'y2', [], 'z2', []);

idx1 = 1;
for i = 1:length(scan)
    for j = 1:length(scan(i).obs)
        all_obs(idx1).orig_idx = idx1;
        all_obs(idx1).scan     = i;
        all_obs(idx1).obs      = j;
        all_obs(idx1).tau      = scan(i).obs(j).obs;
        all_obs(idx1).tau_sig  = scan(i).obs(j).sig;
        all_obs(idx1).i1       = scan(i).obs(j).i1;
        all_obs(idx1).i2       = scan(i).obs(j).i2;
        all_obs(idx1).mjd      = scan(i).mjd;
        all_obs(idx1).e1       = pi/2 - scan(i).stat(scan(i).obs(j).i1).zd;
        all_obs(idx1).e2       = pi/2 - scan(i).stat(scan(i).obs(j).i2).zd;
        all_obs(idx1).az1      = scan(i).stat(scan(i).obs(j).i1).az;
        all_obs(idx1).az2      = scan(i).stat(scan(i).obs(j).i2).az;
        all_obs(idx1).x1       = antenna(scan(i).obs(j).i1).x;
        all_obs(idx1).y1       = antenna(scan(i).obs(j).i1).y;
        all_obs(idx1).z1       = antenna(scan(i).obs(j).i1).z;
        all_obs(idx1).x2       = antenna(scan(i).obs(j).i2).x;
        all_obs(idx1).y2       = antenna(scan(i).obs(j).i2).y;
        all_obs(idx1).z2       = antenna(scan(i).obs(j).i2).z;
        idx1 = idx1 + 1;
    end
end

% Turbulence parameters 
RefFreq    = 8.4e9;             % X-band carrier frequency [Hz]
lambda1Ref = cLight / RefFreq;  % wavelength [m]
Cn21     = 1.2*1e-14; % refractive index structure constant [m^(-2/3)]
HH1      = 2000;      % effective turbulence height [m]
c_aniso1 = 0.01;      % anisotropy parameter c (vertical), a=b=1
a1       = 1.0;       % horizontal anisotropy factor a
u1       = 8.0;       % wind speed [m/s]
u2       = u1;
k_wave1  = 2*pi / lambda1Ref;  % electromagnetic wavenumber [rad/m]
nu1      = 5/6;       % Matern smoothness parameter
L01      = 6000;      % outer scale length [m]
k01      = 2*pi / L01;

NN = totalObservations;
if NN == 0
    error('No observations available.');
end

ManyObs = 0;
if NN > 15000
    ManyObs = 1;
end

% Block-sort observations by station pair (i1, i2)
[~, sortIdx]    = sortrows([[all_obs.i1]', [all_obs.i2]']);
all_obs_blocked = all_obs(sortIdx);
Sigma_turb_var_kov_blocked = zeros(NN, NN);
omega = 2*pi*RefFreq;

if ManyObs == 0

    pref_diag_dynamic = 0.782 * 100 * (k_wave1^2) * HH1 * Cn21 * c_aniso1 * k01^(-5/3);
    k01_1 = 2*pi/L01;
    k01_2 = 2*pi/L01;
    alpha_1 = k01_1 * u1 / a1;
    alpha_2 = k01_2 * u2 / a1;

    % Loop for the calculation of all variances and covariances
    for i = 1:NN
        for j = i:NN

            % Simple mapping function: 1/sin(elevation) 
            s1 = 1 / sin(all_obs_blocked(i).e1); % obs i - station 1
            s2 = 1 / sin(all_obs_blocked(i).e2); % obs i - station 2
            s3 = 1 / sin(all_obs_blocked(j).e1); % obs j - station 1
            s4 = 1 / sin(all_obs_blocked(j).e2); % obs j - station 2

            if i == j
                % Variance
                Sigma_turb_var_kov_blocked(i,i) = (pref_diag_dynamic * (s1^2)) + (pref_diag_dynamic * (s2^2));

            else
                % Time separation [s] between observations i and j
                tau_ij = abs(all_obs_blocked(i).mjd - all_obs_blocked(j).mjd) * 86400;
                if tau_ij == 0, tau_ij = 1e-16; end
            
                valTurb = 0;
                % Obs pairs with time separation > 900 [s] or no shared station get cov = 0
                if tau_ij > 900
                    Sigma_turb_var_kov_blocked(i,j) = valTurb;
                    Sigma_turb_var_kov_blocked(j,i) = valTurb;
                    continue;
                end

                if ManyObs == 1
                    Sigma_turb_var_kov_blocked(i,j) = valTurb;
                    Sigma_turb_var_kov_blocked(j,i) = valTurb;
                    continue;
                end

                sta1 = all_obs_blocked(i).i1;
                sta2 = all_obs_blocked(i).i2;
                sta3 = all_obs_blocked(j).i1;
                sta4 = all_obs_blocked(j).i2;

                dij_1 = 0; dij_2 = 0;

                % Ray separation at turbulence height H (shared stations only)
                if (sta1 == sta3 && sta2 == sta4)
                    dij_1 = dist_rays(all_obs_blocked(i).e1,all_obs_blocked(i).az1,all_obs_blocked(i).x1,all_obs_blocked(i).y1,all_obs_blocked(i).z1,all_obs_blocked(j).x1,all_obs_blocked(j).y1,all_obs_blocked(j).z1,HH1);
                    dij_2 = dist_rays(all_obs_blocked(i).e2,all_obs_blocked(i).az2,all_obs_blocked(i).x2,all_obs_blocked(i).y2,all_obs_blocked(i).z2,all_obs_blocked(j).x2,all_obs_blocked(j).y2,all_obs_blocked(j).z2,HH1);
                elseif (sta1 == sta4 && sta2 == sta3)
                    dij_1 = dist_rays(all_obs_blocked(i).e1,all_obs_blocked(i).az1,all_obs_blocked(i).x1,all_obs_blocked(i).y1,all_obs_blocked(i).z1,all_obs_blocked(j).x2,all_obs_blocked(j).y2,all_obs_blocked(j).z2,HH1);
                    dij_2 = dist_rays(all_obs_blocked(i).e2,all_obs_blocked(i).az2,all_obs_blocked(i).x2,all_obs_blocked(i).y2,all_obs_blocked(i).z2,all_obs_blocked(j).x1,all_obs_blocked(j).y1,all_obs_blocked(j).z1,HH1);
                elseif (sta1 == sta3)
                    dij_1 = dist_rays(all_obs_blocked(i).e1,all_obs_blocked(i).az1,all_obs_blocked(i).x1,all_obs_blocked(i).y1,all_obs_blocked(i).z1,all_obs_blocked(j).x1,all_obs_blocked(j).y1,all_obs_blocked(j).z1,HH1);
                elseif (sta1 == sta4)
                    dij_1 = dist_rays(all_obs_blocked(i).e1,all_obs_blocked(i).az1,all_obs_blocked(i).x1,all_obs_blocked(i).y1,all_obs_blocked(i).z1,all_obs_blocked(j).x2,all_obs_blocked(j).y2,all_obs_blocked(j).z2,HH1);
                elseif (sta2 == sta3)
                    dij_2 = dist_rays(all_obs_blocked(i).e2,all_obs_blocked(i).az2,all_obs_blocked(i).x2,all_obs_blocked(i).y2,all_obs_blocked(i).z2,all_obs_blocked(j).x1,all_obs_blocked(j).y1,all_obs_blocked(j).z1,HH1);
                elseif (sta2 == sta4)
                    dij_2 = dist_rays(all_obs_blocked(i).e2,all_obs_blocked(i).az2,all_obs_blocked(i).x2,all_obs_blocked(i).y2,all_obs_blocked(i).z2,all_obs_blocked(j).x2,all_obs_blocked(j).y2,all_obs_blocked(j).z2,HH1);
                end

                X_1 = alpha_1 * tau_ij;
                X_2 = alpha_2 * tau_ij;
                matern_term_1 = (X_1^(nu1)) * besselk(nu1, X_1);
                matern_term_2 = (X_2^(nu1)) * besselk(nu1, X_2);

                pref_off_1 = 0.7772 * k_wave1^2 * HH1 * Cn21 * c_aniso1 * k01_1^(-5/3);
                pref_off_2 = 0.7772 * k_wave1^2 * HH1 * Cn21 * c_aniso1 * k01_2^(-5/3);

                valTurb = 0;

                % Covariance: only if ray separation < 10000 [m]
                if (sta1 == sta3 && sta2 == sta4)
                    valTurb = (pref_off_1*matern_term_1 * (s1 * s3)) + (pref_off_2*matern_term_2 * (s2 * s4));
                    if dij_1 > 10000 || dij_2 > 10000, valTurb = 0; end
                elseif (sta1 == sta4 && sta2 == sta3)
                    valTurb = (pref_off_1*matern_term_1 * (s1 * s4)) + (pref_off_2*matern_term_2 * (s2 * s3));
                    if dij_1 > 10000 || dij_2 > 10000, valTurb = 0; end
                elseif (sta1 == sta3)
                    valTurb = (pref_off_1*matern_term_1 * (s1 * s3));
                    if dij_1 > 10000, valTurb = 0; end
                elseif (sta1 == sta4)
                    valTurb = (pref_off_1*matern_term_1 * (s1 * s4));
                    if dij_1 > 10000, valTurb = 0; end
                elseif (sta2 == sta3)
                    valTurb = (pref_off_2*matern_term_2 * (s2 * s3));
                    if dij_2 > 10000, valTurb = 0; end
                elseif (sta2 == sta4)
                    valTurb = (pref_off_2*matern_term_2 * (s2 * s4));
                    if dij_2 > 10000, valTurb = 0; end
                else
                    valTurb = 0;
                end

                % Covariance, symmetric fill
                Sigma_turb_var_kov_blocked(i,j) = valTurb;
                Sigma_turb_var_kov_blocked(j,i) = valTurb;
            end
        end
    end

else % use a faster way to calculate only variances

    e1_vals = arrayfun(@(x) x.e1, all_obs_blocked);
    e2_vals = arrayfun(@(x) x.e2, all_obs_blocked);

    % Simple Mapping function
    s1 = 1 ./ sin(e1_vals);
    s2 = 1 ./ sin(e2_vals);

    pref_diag_dynamic = 0.782 * 100 * (k_wave1^2) * HH1 * Cn21 * c_aniso1 * k01^(-5/3);

    % Diagonal Variances
    diag_elements = (pref_diag_dynamic * (s1.^2)) + (pref_diag_dynamic * (s2.^2));
    Sigma_turb_var_kov_blocked = diag(diag_elements);

end

% Back-sort to original observation order 
[~, unsortIdx] = sort(sortIdx);
Sigma_turb_var_kov = Sigma_turb_var_kov_blocked(unsortIdx, unsortIdx); % unit: [rad^2]

addnoise = 1;

% Observations [cm^2] 
Sigma_turb_obs = diag(((mi_observ.^2) .* cLight^2 * 100^2));

% Unit [rad^2] --> [cm^2] 
Sigma_turb_var_kov_cm2 = Sigma_turb_var_kov .* (cLight^2 * 100^2 / omega^2);

% Sum = observations + turbulence matrix
Sigma_turb_all = Sigma_turb_obs + Sigma_turb_var_kov_cm2;


if ManyObs == 0

    % Ensure positive definiteness (with cholesky factorization) and maybe do spectral repair
    [~, p] = chol(Sigma_turb_all);

    if p > 0
        Sigma_turb_all = (Sigma_turb_all + Sigma_turb_all') / 2;
        [V, D] = eig(Sigma_turb_all);
        evals  = diag(D);
        tol    = 1e-9 * max(evals);
        evals(evals < tol) = tol;
        Sigma_turb_all = V * diag(evals) * V';
        Sigma_turb_all = (Sigma_turb_all + Sigma_turb_all') / 2;
        fprintf(' \t Be careful! The variance-covariance matrix has been successfully made positive definite (spectral). \n');
    end

    if p == 0
        fprintf('\t The variance-covariance matrix is positive definite. \n');
    else
        fprintf('\t WARNING: The variance-covariance matrix is NOT positive definite. \n');
    end

else
    fprintf('\t The variance-covariance matrix is positive definite. \n');
    fprintf('\t WARNING: Covariances are NOT calculated due to too many observations (>15000). \n');
    fprintf('\t WARNING: Only variances are calculated! \n');
end

endTime = round(toc);
fprintf('\t Turbulence model calculation: finished after %.1f seconds \n', endTime')

% Ray separation distance at turbulence height H 
function dist1 = dist_rays(el,az,x1,y1,z1,x2,y2,z2,H)

r1 = [x1;y1;z1];
r2 = [x2;y2;z2];

d1 = r1/norm(r1);
d2 = r2/norm(r2);

r1 = r1 + d1*H;
r2 = r2 + d2*H;

bl = r1-r2;

kx = cos(el)*cos(az);
ky = cos(el)*sin(az);
kz = sin(el);

k_hom = [kx;ky;kz;1];

beta = acos(r1(3)/norm(r1));

tmp = r1;
tmp(3)=0;

delta = pi-acos(tmp(1)/norm(tmp));

R1 = [ cos(beta) 0 sin(beta) 0;
            0    1     0     0;
      -sin(beta) 0 cos(beta) 0;
            0    0     0     1 ]; 

R2 = [ cos(delta) -sin(delta) 0 0;
       sin(delta)  cos(delta) 0 0;
            0           0     1 0;
            0           0     0 1];

M = R2' * R1';

ks = M*k_hom;

k = ks(1:3);

tau = dot(bl,k);

dist1 = sqrt(max(norm(bl)^2 - tau^2,0));

end

end