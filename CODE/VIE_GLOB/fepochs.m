% ************************************************************************
%   Description:
%   Find the correct epochs for station coordinates if the velocity was
%   fixed
%
%   Input:										
%       refantbr            structure with information about
%                           discontinuities in station positions
%       refnamec            station names
%       brstart             start of intervals
%       brend               end of intervals
%
%
%   Output:                
%      epoch               epoch for each station interval if the velocity
%                          was fixed --> same as a priori
%
%   External calls: 	
%      -               					    											
%       
%   Coded for VieVS: 
%   20 Dec 2012 by Hana Krásná
%
%   Revision: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

function epoch = fepochs(refantbr,refnamec, brstart, brend)


k=1;
nst=size(refnamec,1);
for i=1:nst
    aname=refnamec(i,:);
    
    estint=[brstart(i) brend(i)];
    for j=1:length(refantbr)
        if refantbr(j).name==aname
            ia=j;
            break
        end
    end

    % find apriori intervals included in the estimation intervals
    aprint=refantbr(ia).break_apr;
    
    % intervals beginn at the same time
    id1 = find(estint(1) <= aprint  & aprint  <= estint(2));
    
    if ~isempty(id1)
        if isempty(find(estint(1) == aprint));
            id1=[id1(1)-1 id1];
        end
        if isempty(find(estint(2) == aprint));
            id1=[id1 id1(end)+1 ];
        end
    end
    
    if isempty(id1)
        xx = find(estint(1) > aprint );
        id1(1)=xx(end);
        clear xx
        xx = find(estint(2) < aprint );
        id1(2)=xx(1);
        
    end
    idAPR = id1(1:end-1);
    
    epoch(i) = refantbr(ia).epoch(idAPR(1));
    % the analysist should check on his own if the epoch yields resonable
    % values, be aware of the epochs in the a priori catalogue!!!
    
end