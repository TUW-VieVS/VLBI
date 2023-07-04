% Return/save to file baseline lengths and proper formal errors!
%
%
% INPUT
%  process_list    process_list (variable or path to mat file). Can be
%                  empty if all files (x_...) are given
%  subdir          LEVEL3 and LEVEL3 subdirectory
%  varargin
%   {1} outfile    full out path of output file. 
%                  If []: no file is output.
%   {2} x_files    struct (length = nSessions) containing x_ files. 
%                  Needed when files are already loaded (eg VieVS)
%                  [] -> not given
%   {3} ant_files  struct (length = nSessions) containing antenna files
%   {4} atpa_files struct (length = nSessions) containing atpa_ files
%   {5} opt_files  struct (length = nSessions) containing apt_ files
%   {6} rigFormErr Rigorous treatment of formal errors. 1=rigorous,
%                  0=simple calculation. Default: 1;
%
% OUTPUT
%  varargout
%   {1} bas      struct, containing baseline information (length = nBas)
%          .ap   a priori baseline length (meters)
%          .corr corrected baseline length (meters)
%          .name baseline name, eg. 'KOKEE   -WESTFORD' (1:8-1:8)
%          .mjd  epoch of baseline (modified julain date)
%          .mbas formal error from least-squares adjstment (rigorously
%                calculated /except: rigFormErr==0) (meters)
%
% CHANGES
% 03.10.2011 K. Teke: Formal errors of the baselines are added
% 30.11.2015 M. Madzak: Added varargout -> Baselines struct is output.
%                       Output can now be textfile or variable or both.
% 11.02.2016 M. Madzak: Name of output file changed (to basel_...)
% 2016-09-21, A. Girdiuk: warning messages are added in case of absence of
%                           station coordinates estimates
% 2017-05-30, A. Girdiuk: opt-file directory info will be used to detect LEVEL1 directory, if the directory with the same as in LEVEL3 does not exist
%                           it will be crased if there is still nothing
% 2018-08-01, D. Landskron: slightly adapted so it can also read vgosDB session names
% 2018-08-14, D. Landskron: some more adaptions for reading vgosDB session names

function [varargout] = bas_out(process_list,subdir,varargin)

%% GET INPUTS
% if a mat file is given instead of already loaded process list
if ~isempty(process_list)
    if strcmpi(process_list(end-3:end), '.mat')
        temp=load(process_list);
        fn=fieldnames(temp);
        process_list=temp.(fn{1});
    end

    pl=size(process_list);
    nSes=pl(1);
else
    nSes=length(varargin{2}{1});
end

if ~isempty(subdir)
	if strcmpi(subdir(end), '/') || strcmpi(subdir(end), '\')
		subdir=subdir(1:end-1);
	end
end

% is textfile written?
out2file=1;

% is variable returend?
out2var=0;
if nargout>0
    out2var=1;
end

if nargin>2
    if isempty(varargin{1})
        out2file=0;
    else
        outfile=varargin{1};
    end
else
	outfile=['../OUT/basel_',subdir,'.txt'];
end

x_filesGiven=0;
antFilesGiven=0;
atpaFilesGiven=0;
optFilesGiven=0;
if nargin>3
    if ~isempty(varargin{2})
        x_files=varargin{2};
        x_filesGiven=1;
    end
    if nargin>4
        if ~isempty(varargin{3})
            ant_files=varargin{3};
            antFilesGiven=1;
        end
        if nargin>5
            if ~isempty(varargin{4})
                atpa_files=varargin{4};
                atpaFilesGiven=1;
            end
            if nargin>6
                if ~isempty(varargin{5})
                    opt_files=varargin{5};
                    optFilesGiven=1;
                end
            end
        end
    end
end

% rigorous treatment of formal errors
rigFormErr=1; % Default:
if nargin>7
    if varargin{6}==0
        rigFormErr=0;
    end
end
    

path= '../';

%% print header to textfile
if out2file==1
    fid = fopen(outfile,'wt');
    %------------------------------------------------
    fprintf(fid,'%%************************************************************\n');
    fprintf(fid,'%% Columns:\n%%\t 1     .... session\n%%\t 2     .... reference time\n%%\t 3     .... baselines\n%%\t 4     .... a priori baseline lengths\n%%\t 5     .... estimated baseline lengths\n%%\t 6     .... formal errors\n'); 
    fprintf(fid,'%% all units are in meters \n');
    fprintf(fid,'%%************************************************************\n');
    fprintf(fid,'%%\n'); 
end

%% preallocate
bas=struct('ap', [], 'corr', [], 'name', [], 'mjd', [], 'mbas', []);
kbas=0;

%% for all sessions
for ip = 1:nSes
    if ~isempty(process_list)
%        sname = process_list(ip,6:end);   % adapted so it can read vgosDB session names as well
        ibckslsh = strfind(process_list(ip,:),'/');
        sname = process_list(ip,ibckslsh(end)+1:end); % adapted for SIMULATED FILES       
        sname_test = strsplit(sname,' ');
        if length(sname_test)>1   % this means that it is a vgosDB session
            sname = sname_test{1};
        end
            
    else
        sname='XXXXXX'; % better but does not work for older data: x_files{1}(ip).session;
    end
   
    if x_filesGiven==1
        x_=x_files{1}(ip);
    else
        load(strcat(path,'DATA/LEVEL3/',subdir,'/x_',num2str(sname)));
    end
    if rigFormErr==1
        if  atpaFilesGiven==1
            atpa_=atpa_files{ip};
        else
            load(strcat(path,'DATA/LEVEL3/',subdir,'/atpa_',num2str(sname)));
        end
    end
    if optFilesGiven==1
        opt_=opt_files{ip};
    else
        load(strcat(path,'DATA/LEVEL3/',subdir,'/opt_',num2str(sname)));
    end
    
	if antFilesGiven
        antenna=ant_files{ip};
    else
        load(strcat(path,'DATA/LEVEL3/',subdir,'/',num2str(sname),'_antenna'));
    end
    
    if isempty([x_.coorx.col]) || isempty([x_.coory.col]) || isempty([x_.coorz.col])
        fprintf('You need to estimate the station coordinates first for session %s!\n',sname);
    else
        index = find(~cellfun(@isempty,{x_.coorx.col}));

            
        index_first = index(1);
        mjd = x_.coorx(index_first).mjd;

        % number of stations
        nstat = length(x_.antenna);

        for k =  1:nstat
            stnam = x_.antenna(k).name;
            % calculate catalogue coordinates at time of interest
            [xa,ya,za] = corap(antenna,stnam,mjd);
            stat.apr(k,:)=[xa,ya,za];

            staest = [x_.coorx(k).val/100,x_.coory(k).val/100,x_.coorz(k).val/100]; %[m]
            if isempty(staest)
                rigFormErr=0; % fixed stations are missing in the N-matrix
                stat.est(k,:) = [0 0 0]; 
            else
                stat.est(k,:) = staest; 
            end

            stat.name(k,:) = stnam;

            stamx=sqrt(((x_.coorx(k).mx/100)^2+(x_.coory(k).mx/100)^2+(x_.coorz(k).mx/100)^2)/3);
            % non-rigourous(!) standard deviation. (in meters)
            if isempty(stamx)
                stat.mx(k,:)=0;
            else
                stat.mx(k,:)=stamx;
            end
        end

        % calculate corrected values as apriori + estimated corrections 
        stat.corr = stat.apr + stat.est;

        % calculate baselines
        % number of baselines
        nbas = nstat*(nstat-1)/2;

        if rigFormErr==1
            % Co-variance matrix
            N = [atpa_.mat];
            Qx = inv(N);
        end

        % calculate baselines and write them to the output file 
        for ib = 1 : nstat-1
           for j = ib+1 : nstat
                stat1 = stat.name(ib,:);
                stat2 = stat.name(j,:);
                % length of apriori baseline [m]
                ap = norm(stat.apr(j,:)-stat.apr(ib,:));
                % length of corrected baseline [m]
                corr = norm(stat.corr(j,:)-stat.corr(ib,:));
                if rigFormErr==1
                    a1 = -(stat.corr(j,1)-stat.corr(ib,1))/corr; 
                    a2 = -(stat.corr(j,2)-stat.corr(ib,2))/corr;
                    a3 = -(stat.corr(j,3)-stat.corr(ib,3))/corr;
                    a4 = -a1;
                    a5 = -a2;
                    a6 = -a3;
                    A = [a1 a2 a3 a4 a5 a6];
                    qx1 = x_.coorx(ib).col;  qy1 = x_.coory(ib).col; qz1 = x_.coorz(ib).col;
                    qx2 = x_.coorx(j).col;   qy2 = x_.coory(j).col;  qz2 = x_.coorz(j).col;
                    Qd = [];
                    Qd = [Qx(qx1,qx1) Qx(qx1,qy1) Qx(qx1,qz1) Qx(qx1,qx2) Qx(qx1,qy2) Qx(qx1,qz2)
                          Qx(qy1,qx1) Qx(qy1,qy1) Qx(qy1,qz1) Qx(qy1,qx2) Qx(qy1,qy2) Qx(qy1,qz2)
                          Qx(qz1,qx1) Qx(qz1,qy1) Qx(qz1,qz1) Qx(qz1,qx2) Qx(qz1,qy2) Qx(qz1,qz2)
                          Qx(qx2,qx1) Qx(qx2,qy1) Qx(qx2,qz1) Qx(qx2,qx2) Qx(qx2,qy2) Qx(qx2,qz2)
                          Qx(qy2,qx1) Qx(qy2,qy1) Qx(qy2,qz1) Qx(qy2,qx2) Qx(qy2,qy2) Qx(qy2,qz2)
                          Qx(qz2,qx1) Qx(qz2,qy1) Qx(qz2,qz1) Qx(qz2,qx2) Qx(qz2,qy2) Qx(qz2,qz2)];  
                    md = sqrt(A*Qd*A')*opt_.mo*0.01; % [m] 
                else
                    md=sqrt( (stat.mx(ib))^2+(stat.mx(j))^2 ); % meters
                end
                if out2file==1
                    fprintf(fid,'%s %s %s - %s %15.4f %15.4f %7.4f\n',num2str(sname),num2str(mjd),stat1(1:8),stat2(1:8),ap,corr,md);
                end

                % output to variable
                if out2var==1
                    curBasN=[stat1,'-',stat2];
                    curBasN2=[stat2,'-',stat1];

                    % see if we already have that baseline
                    foundBasLog=strcmpi({bas.name},curBasN) | ...
                        strcmpi({bas.name},curBasN2);
                    if sum(foundBasLog)==0 % we have a new baseline
                        kbas=kbas+1;
                        toBasInd=kbas;
                        toBasFieldsInd=1;
                        bas(toBasInd).name=curBasN;
                    else
                        toBasInd=find(foundBasLog);
                        toBasFieldsInd=length(bas(toBasInd).mjd)+1;
                    end
                    bas(toBasInd).ap(toBasFieldsInd)=ap;
                    bas(toBasInd).corr(toBasFieldsInd)=corr;
                    bas(toBasInd).mjd(toBasFieldsInd)=mjd;
                    bas(toBasInd).mbas(toBasFieldsInd)=md;
                end
            end
        end
       
    end
end

if out2file==1
    fclose(fid);
end

% assign output: variable 'bas'
if out2var==1
    varargout{1}=bas;
end

