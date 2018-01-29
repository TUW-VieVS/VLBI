function [itim,idoy]=dday(iy,im,id,ih,imi)

idoy=zeros(size(iy));
itim=idoy;

for ik=1:length(iy)
    daymon=[0 31 28+(mod(iy(ik),4)==0) 31 30 31 30 31 31 30 31 30 31];
    idoy(ik)=sum(daymon(1:im(ik)))+id(ik);
    itim(ik)=round((ih(ik)+imi(ik)/60)*10);
end