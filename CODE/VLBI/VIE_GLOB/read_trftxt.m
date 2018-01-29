% ************************************************************************
%   Description:
%   This function reads a TRF catalogue in .txt format
%
%   Input:										
%       trffile   path to the .txt file
%                           
%   Output:                
%       trf       trf catalogue saved in a matlab structure
%
%   External calls: 	
%      -               					    											
%       
%   Coded for VieVS: 
%    Hana Spicakova: the code was copied from vie_init/read_ngs.m prepared by T.Nilsson
%
%%

function trf=read_trftxt(trffile)

fid=fopen(trffile,'r');
nt=0;
trf=[];
while ~feof(fid)
    str=fgetl(fid);
    if length(str)>10
        if (str(1)~='%')
            nuffs=str2num(str(9:length(str)));
            if length(nuffs)>=9
                numt=0;
                for a=1:length(trf)
                    if strcmp(trf(a).ivsname,str(1:8))
                        numt=a; break;
                    end
                end
                if numt==0
                    nt=nt+1;
                    numt=nt;
                    trf(numt).ivsname=str(1:8);
                    trf(numt).break=[];
                end
            else
                nt=nt+1;
                numt=nt;
                trf(numt).ivsname=str(1:8);
                trf(numt).break=[];
            end
            nbr=length(trf(numt).break)+1;

            trf(numt).break(nbr).x=nuffs(1);
            trf(numt).break(nbr).y=nuffs(2);
            trf(numt).break(nbr).z=nuffs(3);
            if length(nuffs)>=9
                trf(numt).break(nbr).start=nuffs(8);
                trf(numt).break(nbr).end=nuffs(9);
            else
                trf(numt).break(nbr).start=0;
                trf(numt).break(nbr).end=999999;
            end
            if length(nuffs)>=7
                trf(numt).break(nbr).epoch=nuffs(7);
                trf(numt).break(nbr).vx=nuffs(4);
                trf(numt).break(nbr).vy=nuffs(5);
                trf(numt).break(nbr).vz=nuffs(6);
            else
                trf(numt).break(nbr).epoch=51544;
                trf(numt).break(nbr).vx=0;
                trf(numt).break(nbr).vy=0;
                trf(numt).break(nbr).vz=0;
            end
        end
    end
end
fclose(fid);

