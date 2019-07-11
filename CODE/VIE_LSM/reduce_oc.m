% ************************************************************************
%   Description:
%   function to initial treatment for the clock erros before the main estimation,
%   clock breaks are the used as the boundaries of different clock functions (these
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
%       'fname'
%       'dirpth'
%
%   Output:
%       'oc_observ'         (1,nobserv)         o-c vector of real observations
%       'first_solution'    structure array     several outputs of first solution for book keeping
%       'opt'               structure array
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
%       breaks.
%   14 Jun 2012 by Matthias Madzak: Saving first solution residuals to
%       file in LEVEL3. Therefore the user-subdir and the sessionname (fname)
%       are given as parameters.
%   26 Apr 2013 by Sigrid Boehm: screen output of m0 changed to chi^2
%   10 Oct 2016 by A. Hellerschmied: fields "obs_type" and "allSatelliteNames" added to "res" structure
%   01 Mar 2017 by A. Hellerschmied: - "Manually find clock breaks" function removed
%                                    - Code revised
% ************************************************************************
function [oc_observ,first_solution,opt] = reduce_oc(nobserv,na,ntim,scan,Pobserv,oc_observ,opt,per_stat,fname,dirpth)

numberOfLSMs = size(oc_observ,2);


% Forming clock breaks vector
breaks = [];
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
a_offset    = []; 
a_rate      = []; 
a_quad      = [];
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
first_solution.mo = sqrt(diag((first_solution.v_clk'*Pobserv*first_solution.v_clk)) / (size(A,1)-length(first_solution.clk_val))); % [cm]
disp(['chi-squared of first solution: ',num2str(first_solution.mo'.^2)]);
nP = length(Q_clk);
first_solution.clk_mx = repmat(first_solution.mo',nP,1).*repmat(sqrt(diag(Q_clk)),1,numberOfLSMs); % std. dev. of estimated clock offsets [cm]

% Residuals of first solution
oc_observ_first = oc_observ - A * [first_solution.clk_val];

% Reduction of clock functions (first treatment); the zenith delays are kept
oc_observ = oc_observ - A(:,1:colclk) * first_solution.clk_val(1:colclk,:);


% +++ save first-solution residuals to variable +++
% Source names:
res.allStatNames  = {opt.stat.name}';
if isfield(opt, 'source')
    if ~isempty(opt.source)
        res.allSourceNames = {opt.source.name}';
    else
        res.allSourceNames = {};
    end
else
    res.allSourceNames = {};
end
if isfield(opt, 'satellite')
    if ~isempty(opt.satellite)
        res.allSatelliteNames = {opt.satellite.name}';
    else
        res.allSatelliteNames = {};
    end
    res.allSatelliteNames = {opt.satellite.name}';
else
    res.allSatelliteNames = {};
end

% Init:
res.baselineOfObs   = zeros(nobserv,2);
res.source          = zeros(nobserv,1);
res.mjd             = zeros(nobserv,1);

curInd              = 1;
for iScan = 1 : size(scan, 2)
    % get station indices of current scan
    statOfCurScan = [[scan(iScan).obs.i1]', [scan(iScan).obs.i2]'];
    nObsOfThisScan = size(statOfCurScan,1);
    % source index of current scan
    res.source(curInd:nObsOfThisScan+curInd-1) = repmat(scan(iScan).iso, nObsOfThisScan,1);
    % obs type
    res.obs_type(curInd:nObsOfThisScan+curInd-1) = repmat(scan(iScan).obs_type, nObsOfThisScan,1);
    % get baselines, mjd and sources of observations
    res.baselineOfObs(curInd:nObsOfThisScan+curInd-1,1:2) = statOfCurScan;
    res.mjd(curInd:nObsOfThisScan+curInd-1) = repmat(scan(iScan).mjd, nObsOfThisScan,1);
    curInd = curInd + nObsOfThisScan;
end
res.obs_type = res.obs_type';
res.firstVal = oc_observ_first * -1;

% the following is correct (same as in ngs file!)

% save
if ~isempty(dirpth) && ~exist(['../DATA/LEVEL3/',dirpth], 'dir')
    mkdir(['../DATA/LEVEL3/',dirpth])
end

if numberOfLSMs == 1
    save(['../DATA/LEVEL3/', dirpth, '/res_', fname, '.mat'], 'res');
end

