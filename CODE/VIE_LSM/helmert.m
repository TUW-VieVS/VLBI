% ************************************************************************
%   Description:
%   function to form NNT/NNR and NNS condition equations
%   for datum definition of TRF
%
%   Reference:
%       Hana Krasna, thesis: "Estimation of solid Earth tidal parameters
%                              and FCN with VLBI"
%
%   Input:
%       'n_'         structure array            number of estimates (pwlo or one offset)
%       'na'         (1,1)                      number of antennas
%       'xo'         (1,na)                     apriori TRF (x) coordinates of all antennas in the session
%       'yo'         (1,na)                     apriori TRF (y) coordinates of all antennas in the session
%       'zo'         (1,na)                     apriori TRF (z) coordinates of all antennas in the session
%       'opt'        structure array            (for info. /DOC/opt.doc)
%       'sum_dj'     (1,number of models)       total number of estimates for each included model
%       'N'          (sum_dj(end),sum_dj(end))  datum free normal equation matrix
%
%   Output:
%       'N' (sum_dj(end)+cond,sum_dj(end)+cond) normal equation matrix with NNT/NNR condition equations
%
%   External calls:
%
%   Coded for VieVS:
%   12 May 2009 by Kamil Teke
%
%   Revision:
%   06 Dec 2009 by Kamil Teke: header added
%   30 Nov 2016 by A. Girdiuk: function reviewed, non-constrained solution allowed
%   08 May 2017 by A. Hellerschmied: Changes for estimating satellite postion offsets (pwl)
% ************************************************************************
function [N] = helmert(n_,na,xo,yo,zo,opt,sum_dj,N) % Generalized Inverse for antenna coordinates

B1 = zeros(7,na*3);
cc = 1/sqrt(xo*xo'+yo*yo'+zo*zo');
xii = cc*xo; yii = cc*yo; zii = cc*zo;

len_level=0;

for istat = 1 : na
    %--------------------
    B_istat = [  1           0             0
        0           1             0
        0           0             1
        0       -zii(istat)  yii(istat)
        zii(istat)   0         -xii(istat)
        -yii(istat)  xii(istat)       0
        xii(istat)  yii(istat)  zii(istat)
        ];
    
    for iter = 0 : n_(istat).xyz-1
        B1( : , istat + iter + len_level)                   = B_istat(:,1); % X
        B1( : , istat + iter + sum([n_.xyz])   + len_level) = B_istat(:,2); % Y
        B1( : , istat + iter + sum([n_.xyz])*2 + len_level) = B_istat(:,3); % Z
    end
    len_level= len_level + n_(istat).xyz-1;
end

nnt = [opt.stat.nnt_inc]==0;
nnr = [opt.stat.nnr_inc]==0;
nns = [opt.stat.nns_inc]==0;

fprintf('Station not in NNT: ')
for i=1:na
    if nnt(i)
        fprintf('%s ',opt.stat(i).name)
    end
end
fprintf('\n')
fprintf('Station not in NNR: ')
for i=1:na
    if nnr(i)
        fprintf('%s ',opt.stat(i).name)
    end
end
fprintf('\n')
fprintf('Station not in NNS: ')
for i=1:na
    if nns(i)
        fprintf('%s ',opt.stat(i).name)
    end
end
fprintf('\n')

if ~any(nns==0)
    k=6;
    B1(7,:)=[];
else
    k=7;
    
    B1(7,nns)=0;
    B1(7,find(nns)+na)=0;
    B1(7,find(nns)+na*2)=0;
end

B1(:,nnt)                =zeros(1,k,sum(nnt));
B1(:,find(nnt)+na)       =zeros(1,k,sum(nnt));
B1(:,find(nnt)+na*2)     =zeros(1,k,sum(nnt));

B1(:,nnr)                =zeros(1,k,sum(nnr));
B1(:,find(nnr)+na)       =zeros(1,k,sum(nnr));
B1(:,find(nnr)+na*2)     =zeros(1,k,sum(nnr));


B1(~any(B1,2),:) = [];

if sum_dj(13) ~= 0
    if size(B1,1)~=0
        B2 = zeros(size(B1,1),sum_dj(13));
        if length(sum_dj) >= 16     % satellite positions pwl
            B  = horzcat(B2,B1);
            B  = horzcat(B, zeros(size(B, 1), (size(N, 2) - size(B, 2))));
        else
            B  = horzcat(B2,B1);
        end
    else
        B = B1;
    end
else
    B = B1;
end

if ~isempty(B)
    K = zeros(size(B,1),size(B,1));
    N = horzcat(vertcat(N,B),vertcat(B',K));
    fprintf('!!! NNT and NNR conditions are introduced to matrix N for station coordinates!!!\n');
else
    fprintf('Attention! Stations are not constrainted\n');
end

