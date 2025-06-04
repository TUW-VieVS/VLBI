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
%   22 Jan 2025 by H. Wolf: improved code for performance
% ************************************************************************
function [N] = helmert(n_,na,xo,yo,zo,opt,sum_dj,N) % Generalized Inverse for antenna coordinates

    if opt.stc_all == 1
        pos_x = 13; 
    elseif opt.stc_sat == 1
        pos_x = 27; 
    elseif opt.stc_qu == 1
        pos_x = 30;
    elseif opt.stc_qs == 1
        pos_x = 30; % first qu
    end
    
    B1 = zeros(7,na*3);
    cc = 1/sqrt(xo*xo'+yo*yo'+zo*zo');
    xii = cc*xo; yii = cc*yo; zii = cc*zo;
    
    len_level=0;
    names = strings(1, na);
    for istat = 1 : na
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
        names(istat) = opt.stat(istat).name;
    end
    
    nnt = [opt.stat.nnt_inc]==0;
    nnr = [opt.stat.nnr_inc]==0;
    nns = [opt.stat.nns_inc]==0;
    
    not_nnt_names = names(nnt==1);
    not_nnr_names = names(nnr==1);
    not_nns_names = names(nns==1);
    fprintf('Station not in NNT: ')
    fprintf('%s ', not_nnt_names(1:end))
    fprintf('\nStation not in NNR: ')
    fprintf('%s ', not_nnr_names(1:end))
    fprintf('\nStation not in NNS: ')
    fprintf('%s ', not_nns_names(1:end))
    fprintf('\n')
      
    %remove nnt for stations which are not in datum
    B1(1:3,nnt)                =zeros(1,3,sum(nnt));
    B1(1:3,find(nnt)+na)       =zeros(1,3,sum(nnt));
    B1(1:3,find(nnt)+na*2)     =zeros(1,3,sum(nnt));
    
    %remove nnr for stations which are not in datum
    B1(4:6,nnr)                =zeros(1,3,sum(nnr));
    B1(4:6,find(nnr)+na)       =zeros(1,3,sum(nnr));
    B1(4:6,find(nnr)+na*2)     =zeros(1,3,sum(nnr));

    %remove nns for stations which are not in datum
    if ~any(nns==0)
        B1(7,:)=[];
    else   
        B1(7,nns)=0;
        B1(7,find(nns)+na)=0;
        B1(7,find(nns)+na*2)=0;
    end
    
    B1(~any(B1,2),:) = [];
    
    %add zeros before and behind
    if opt.stc_qs==0
        B2 = zeros(size(B1,1),sum_dj(pos_x));
        B  = horzcat(B2,B1);
        B  = horzcat(B, zeros(size(B, 1), (size(N, 2) - size(B, 2))));
    else
        B1_nnr = B1;
        B1_nnr (1:3, :) = zeros(3, size(B1,2));
        B_qs  = horzcat(B1_nnr,B1);
        B2 = zeros(size(B_qs,1),sum_dj(27));
        B  = horzcat(B2,B_qs);
        B  = horzcat(B, zeros(size(B, 1), (size(N, 2) - size(B, 2))));
    end
    
    if ~isempty(B)
        K = zeros(size(B,1),size(B,1));
        N = horzcat(vertcat(N,B),vertcat(B',K));
    end
end