% ************************************************************************
%   Description:
%   In this function the datum conditions are applied and the least squares 
%   adjustment is carried out.
%
%   Input:										
%      A            design matrix
%      Pobserv      weight matrix
%      opt          settings
%      oc_observ    o-c vector of observations
%      n_observ     number of observations
%      sum_dj       number of estimated parameters
%      xo           a-priori coordinate x
%      yo           a-priori coordinate y
%      zo           a-priori coordinate z 
%      n_			number of estimates (pwlo or one offset)
%      na           number of antennas
%      ns_q         number of sources (quasars)
%      ra           right ascension of sources
%      de           declination of sources
%      nistat		the number of reference clock (in the order of "antenna")
%      antenna      antenna struct
%
% 
%   Output:
%     x                 estimates
%     v                 residuals
%     v_real			residuals real
%     Qxx               Qxx matrix
%     N					normal equation matrix
% 
%   External calls: 	
%       helmert.m   
%       
%   Coded for VieVS: 
%   25 Jan 2025 by Helene Wolf: as external function from vie_lsm
%
%   Revision: 
%
%
% ************************************************************************
function [x, v, v_real, Qxx, N] = process(A, Pobserv, opt, oc_observ, n_observ, sum_dj, xo, yo, zo, n_, na, ns_q, ra, de, nistat, antenna)

    N = A'*Pobserv*A;
    N = full(N);
    % apply NNT and/or NNR for station coordinates
    if opt.stc == 1 && opt.addDatumCd == 1
        if opt.nnt_stc == 1 || opt.nnr_stc == 1 || opt.nns_stc == 1
            [N] = helmert(n_, na, xo, yo, zo, opt, sum_dj, N);
        end
    end

    % apply NNR for source coordinates
    if opt.est_sourceNNR==1
        if sum([opt.source.nnr_inc]) > 2
            fprintf('NNR condition for source coordinates is introduced to N matrix!\n');
           % fprintf(1,'Not enough sources for NNR condition. All sources are used for NNR instead.\n');
           % [opt.source.nnr_inc] = deal(1);
           %end     
           [N] = nnr_cond(ns_q,ra,de,opt,sum_dj,N);
        else
            fprintf('Less than 3 datum sources!\n');
       end
    end

    fprintf('clock %s is selected as the ref. clock for the main solution\n',antenna(nistat).name);

    % Condition number of N matrix 
    % [unitless], close to 1 -> well-conditioned matrix
    condn = condest(N);
    [opt.condn] = condn;

    Qxx = inv(N); % [cm2] & [mas2]
    Qxx = Qxx(1:sum_dj(length(sum_dj)),1:sum_dj(length(sum_dj)));

    n = A'*Pobserv*oc_observ; % [1/cm] & [1/mas]
    x = Qxx*n; % [cm] & [mas]

    v = A*x - oc_observ; % [cm]
    v_real = v(1:n_observ,:);
end