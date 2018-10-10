% ************************************************************************
%   Description:
%   Creates B_vties - matrix with velocity ties
%
%
%   Input:
%      refnamec          names of antennas, also multiple if there is a
%                        break in the position
%      velties           structure with names of antennas, where the
%                        velocity should be equal
%      velconst          names of antennas with breals, where only one
%                        velocity will be estimated for all intervals
%   Output:                
%      B_vties           matrix with velocity constraints
%
%
%   Coded for VieVS: 
%   10 Jul 2011 by Hana Spicakova
%
%   Revision: 
%   25 Mar 2013 by Hana Krasna: bug corrected, which appeared if a station
%       with breaks was constrained to station without breaks. (Station with
%       breaks has to be constrained to a constant velocity at all
%       intervals.)
%   21 Jun 2013 by Hana Krasna: it is possible to distinquish between
%              intervals for velocity constraints at one station with
%              breaks
%%


function B_vties = create_B_vties(refnamec,velties,velconst,velconst_interv)

s1=size(refnamec,1);
B_vtie(1,1:s1)=zeros(s1,1);

k=0;
for i=1:length(velties)   
   IndRefnamec=[];
   for j=1:size(velties(i).aname,1)
        f=find(strcmp(cellstr(velties(i).aname(j,:)),cellstr(refnamec))==1);
        if ~isempty(f)
            IndRefnamec = [IndRefnamec f(1)];
        end
   end
   
   for j=2:length(IndRefnamec)
       k=k+1;
       B_vtie(k,IndRefnamec(1))= 1;
       B_vtie(k,IndRefnamec(j))=-1;
   end
end
   
%% constant velocity for stations with breaks

for i = 1:size(velconst,1)
    clear f vn ind1
    f=find(strcmp(cellstr(velconst(i,1:8)),cellstr(refnamec)));
    vn=str2num(velconst_interv{i}); % intervals where the velocity should be constant
    if ~isempty(f)
        if nonzeros(vn)
            ind1=f(vn);
        else
            ind1=f;
        end
    else
        ind1=0;
    end
    
    ind2=nonzeros(ind1);
    if length(ind2)>1
       for ij=2:length(ind2);
           k=k+1;
           B_vtie(k,ind2(1))= 1;
           B_vtie(k,ind2(ij))=-1;
       end
   end
 end



%%
s2=size(B_vtie,1);
B_vties=[B_vtie       zeros(s2,s1) zeros(s2,s1)
         zeros(s2,s1) B_vtie       zeros(s2,s1)
         zeros(s2,s1) zeros(s2,s1) B_vtie ];

     
     
if k==0
   B_vties=[];
end
    
