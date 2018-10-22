% ************************************************************************
%   Description:
%   function to initial treatment for the clock erros before the main estimation,
%   clock breaks are used as the boundaries of different clock functions (these
%   clock functions include only one offset estimate, or one offset + rate,
%   or one offset + one rate + one quadratic term, according to the options)
%   'breaks(istat).clk' structure array specifies the epoch exactly when the
%   clock break occured [firstobstime - clockbreak1 - lastobstime].
%   The estimated clock errors (considering the clock breaks by partitioning
%   the clock function) are reduced from o-c (named as : FIRST CLOCK TREATMENT)
%   Zenith wet delays estimated and reduced from o-c only for ploting
%   pruposes in order to check and see clearly if any missing clock break occured.
%
%   Reference:
%
%   Input:
%       'nobserv'           (1,1)               number observations in a session
%       'na'                (1,1)               number of antennas
%       'ntim'              (1,1)               number of scans in a session
%       'scan'              structure array     (for info. /DOC/scan.doc)
%       'Pobserv'           (nobserv,nobserv)   weights of real observtaions
%       'oc_observ'         (1,nobserv)         o-c vector of real observations
%       'opt'               structure array     (for info. /DOC/opt.doc)
%       'per_stat'          structure array     (for info. /DOC/per_stat.doc)
%
%   Output:
%       'oc_observ'         (1,nobserv)         o-c vector of real observations
%       'first_solution'    structure array     several outputs of first solution for book keeping
%
%   External calls:
%
%   Coded for VieVS:
%   06 July 2009 by Johannes Boehm & Kamil Teke
%
%   Revision:
%   15 Aug 2009 by Kamil Teke: parameters of clock function and formal
%                              errors of these estimates attached to output : 'x_'
%   21 Aug 2009 by Johannes Boehm : zenith delay offsets are also estimated in the first solution,
%                                   but they are not used to reduce the observations
%   06 Dec 2009 by Kamil Teke: header added
%   15 Sep 2010 by Johannes Boehm: The mapping function for the estimation of the zenith delays
%                                  was changed to 1/cos(zd)
%   18 Apr 2011 by Tobias Nilsson: Implemented a tool for detecting clock
%   breaks.
%   23 Feb 2017 by Hana Krasna: param. addnoise added as input
%   01 Mar 2017 by A. Hellerschmied: - "Manually find clock breaks" function removed
%                                    - Code revised
% ************************************************************************
function [oc_observ,first_solution,opt] = reduce_oc_scan(nobserv,na,ntim,scan,oc_observ,opt,c,addnoise)

% ASSIGNING THE INPUT PARAMETERS TO EACH STATION
% creation of per_stat structure
for istat = 1:na % number of stations
    itim    = 0; 
    i       = 0; 
    k       = 0;
    for itim = 1:ntim % number of scans per session
        for iobs = 1:scan(itim).nobs % number of observations per scan
            i = i + 1; % i : i. observation in the session
            i1 = scan(itim).obs(iobs).i1;
            i2 = scan(itim).obs(iobs).i2;
            if i1 == istat || i2 == istat
                k = k + 1; % k : k. observation of the specific station
                per_stat(istat).mjd(k) = scan(itim).mjd; % The times of scans [day]
                per_stat(istat).oc_nob(k) = i;
                per_stat(istat).zd(k) = scan(itim).stat(istat).zd;  % Boehm 21 Aug 2009, 15 Sep 2010
            end
            if i1 == istat
                per_stat(istat).first(k) = -1;
                per_stat(istat).other(k) = i2;  % Boehm 21 Aug 2009
            end
            if i2 == istat
                per_stat(istat).first(k) = +1;
                per_stat(istat).other(k) = i1;  % Boehm 21 Aug 2009
            end
        end
    end
end


% FORMING THE WEIGHT MATRIX OF THE OBSERVATIONS (Pobserv)
temp = [scan.obs];
mi_observ = [temp.sig]; % [seconds]
nmi_observ = sqrt((addnoise/c)^2.+(mi_observ).^2); % [seconds]
so = sqrt(((nmi_observ.*c*100)*(nmi_observ.*c*100)')/nobserv); % [cm] apriori std. dev. of unit weight
opt.so = so;
Pobserv = diag(sparse(1./((nmi_observ.^2).*c^2*100^2))); % [1/cm^2]
Qll = diag(sparse((nmi_observ.^2).*c^2*100^2)); % [cm^2]
fprintf('apriori std. dev. of unit weight. : %3.4f\n',so);


% Forming clock breaks vector
count = 0; breaks = [];
for istat = 1 : na
    if opt.treat_breaks == 1
        if ~isempty(opt.stat(istat).clkbreak)
            breaks(istat).clk = vertcat(scan(1).mjd,opt.stat(istat).clkbreak,scan(ntim).mjd);
        else
            breaks(istat).clk = vertcat(scan(1).mjd,scan(ntim).mjd);
        end
    else
        breaks(istat).clk = vertcat(scan(1).mjd,scan(ntim).mjd);
    end
end

% Forming clock functions according to the clock breaks vector above
a_offset = []; a_rate = []; a_quad = [];
for istat = 1 : na % stations
    for part = 1 : length(breaks(istat).clk)-1 % clock breaks
        a(istat).clk(part).offset(nobserv,1) = 0;
        a(istat).clk(part).rate(nobserv,1) = 0;
        a(istat).clk(part).quad(nobserv,1) = 0;
        for k = 1 : length(per_stat(istat).mjd)
            mjd =  per_stat(istat).mjd(k);
            first = per_stat(istat).first(k);
            if (mjd >= breaks(istat).clk(part)) && (mjd <= breaks(istat).clk(part + 1))
                if opt.firstclock == 0 || opt.firstclock == 1 || opt.firstclock == 2
                    a(istat).clk(part).offset(per_stat(istat).oc_nob(k),1) = first;
                end
                if opt.firstclock == 1 || opt.firstclock == 2
                    a(istat).clk(part).rate(per_stat(istat).oc_nob(k),1) = first * (mjd - breaks(istat).clk(part));
                end
                if opt.firstclock == 2
                    a(istat).clk(part).quad(per_stat(istat).oc_nob(k),1) = first * (mjd - breaks(istat).clk(part))^2;
                end
            end
        end
        if istat ~= opt.ref_first_clk % Not taking reference clock
            if opt.firstclock == 0 || opt.firstclock == 1 || opt.firstclock == 2
                a_offset = horzcat(a_offset,a(istat).clk(part).offset);
            end
            if opt.firstclock == 1 || opt.firstclock == 2
                a_rate = horzcat(a_rate,a(istat).clk(part).rate);
            end
            if opt.firstclock == 2
                a_quad = horzcat(a_quad,a(istat).clk(part).quad);
            end
        end
    end
end

% Forming on zwd offset per station
a_zwd = [];
for istat = 1 : na % stations
    for k = 1 : length(per_stat(istat).mjd)
        first = per_stat(istat).first(k);
        mf = 1.0/cos(per_stat(istat).zd(k));  % Boehm, 15 Sep 2010
        oc_nob = per_stat(istat).oc_nob(k);
        a_zwd(oc_nob,istat) = first*mf;       % Boehm, 15 Sep 2010
    end
end

colclk = size([a_offset,a_rate,a_quad],2); % number of columns for the clocks in the first solution

A = horzcat(a_offset,a_rate,a_quad,a_zwd);

Q_clk = inv(A'*Pobserv*A);
b_clk = A'*Pobserv*oc_observ;

first_solution.ref_st_name = opt.stat(opt.ref_first_clk).name;
first_solution.ref_st = opt.ref_first_clk;
disp(['clock ',first_solution.ref_st_name,' is selected as the ref. clock for the first solution']);
first_solution.breaks = breaks;
first_solution.clk_val = Q_clk*b_clk; % [cm]
first_solution.v_clk = A*first_solution.clk_val - oc_observ; % [cm]
first_solution.mo = sqrt((first_solution.v_clk'*Pobserv*first_solution.v_clk) / (size(A,1)-length(first_solution.clk_val))); % [cm]
disp(['chi-squared of first solution: ',num2str(first_solution.mo)]);
first_solution.clk_mx = first_solution.mo*sqrt(diag(Q_clk)); % std. dev. of estimated clock offsets [cm]

% Residuals of first solution
% oc_observ_plot = oc_observ - A*[first_solution.clk_val];

% Reduction of clock functions (first treatment); the zenith delays are kept
oc_observ = sparse(oc_observ - A(:,1:colclk)*[first_solution.clk_val(1:colclk)]);

