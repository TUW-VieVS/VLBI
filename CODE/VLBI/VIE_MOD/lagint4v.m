% ************************************************************************
%   Description:
%   This subroutine performs lagrangian interpolation within a set of (X,Y)
%   pairs to give the Y value corresponding to XINT. This program uses a 
%   window of 4 data points to perform the interpolation. If the window
%   size needs to be changed, this can be done by changing the indices in
%   the loops for variables m and j. Translated occam subroutine lagint4.f,
%   distributed by Daniel Gambis. 
%
%   Input:										
%      X                   Array of values of the indepnedent variable
%      Y                   Array of Function  values corresponding to X
%      XINT                The X-value for shich estimate of Y is desired
%                
%   Output:
%      YOUT                The Y-value returned
% 
%   External calls: 	
%      cart2phigd.m, ren2xyz.m   
%
%   Coded for VieVS: 
%   13 Jan 2009 by Hana Spicakova
%
%   Revision: 
%   20 Jan 2009 by Lucia Plank, allow vector input for lagint4, new name
%   19 Jul 2011 by Lucia Plank: replace for with find
%   16 Sep 2011 by Hana Spicakova: fixed the bug, going back from +6 days
%   to +5 days
%   12 Mar 2015 by Matthias Madzak: Change of epoch (index) selection. Now
%   always 5 values before and after ar taken (if epoch is equal, value is
%   the same anyway!).
% *************************************************************************


function [YOUTs]=lagint4v(X,Y,XINT)

n=length(XINT);
YOUTs=zeros(n,1);

for s=1:n

    N=length(X);
    YOUT = 0.0; 

%     for i=1:(N-1)
%         if (XINT(s) >= X(i)) && (XINT(s) < X(i+1))
%             k=i;
%         end
%     end
%     ind=find(XINT(s)>=X);
%     k=ind(end);
    
    % new selection of proper intervals (five before , five after -> ten
    % values always)
    indSmallerFlip=flipud(find(X<=XINT(s)));
    indFirst=indSmallerFlip(5);
    indLarger=find(X>=XINT(s));
    indLast=indLarger(5);
    
%     if k < 2
%         k=2;
%     end
%     
%    
%     if k > (N-2)
%         k=N-2;
%     end
    
    for m = indFirst:indLast
         TERM = Y(m);
        for j = indFirst:indLast
            if ( m ~= j) 
                 TERM = TERM * (XINT(s) - X(j))/(X(m) - X(j));
            end 
        end
            YOUT = YOUT + TERM;
    end
    YOUTs(s)=YOUT;
end