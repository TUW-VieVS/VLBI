% Purpose  
%   Analyze the 4-source subnet. 
% History  
%   2010-04-06   Jing SUN   Created
%   


function [catpair] = isub4grid(cat)

% regular tetrahedron
cor1 = [0; 0; 1];
cor2 = [0;         -2*sqrt(2)/3;    -1/3];
cor3 = [ sqrt(6)/3;   sqrt(2)/3;    -1/3];
cor4 = [-sqrt(6)/3;   sqrt(2)/3;    -1/3];
%
gridnum = length(cat);
for igrid1 = 1 : gridnum
    % rangeid1
    rangeid1 = igrid1;
    % ra1, de1
    ra1 = cat(rangeid1).ramean;
    de1 = cat(rangeid1).demean;
    fai1 = ra1;
    sita1 = pi/2 - de1; 
    % rotation matrix
    ry = [cos(-sita1)  0 -sin(-sita1);
          0            1        0;
          sin(-sita1)  0 cos(-sita1)];
    rz = [cos(-fai1) sin(-fai1) 0;
         -sin(-fai1) cos(-fai1) 0;
          0         0         1];
    % ra2, de2
    xyz2 = rz * ry * cor2;
    x2 = xyz2(1);
    y2 = xyz2(2);
    z2 = xyz2(3);
    r2 = sqrt(x2*x2+y2*y2+z2*z2);
    sita2 = acos(z2/r2);
    fai2  = atan2(y2,x2);
    if (fai2 < 0)
        fai2 = fai2 + 2*pi;
    end
    ra2 = fai2;
    de2 = pi/2 - sita2;
    % ra3,de3
    xyz3 = rz * ry * cor3;
    x3 = xyz3(1);
    y3 = xyz3(2);
    z3 = xyz3(3);
    r3 = sqrt(x3*x3+y3*y3+z3*z3);
    sita3 = acos(z3/r3);
    fai3  = atan2(y3,x3);
    if (fai3 < 0)
        fai3 = fai3 + 2*pi;
    end
    ra3 = fai3;
    de3 = pi/2 - sita3;
    % ra4,de4
    xyz4 = rz * ry * cor4;
    x4 = xyz4(1);
    y4 = xyz4(2);
    z4 = xyz4(3);
    r4 = sqrt(x4*x4+y4*y4+z4*z4);
    sita4 = acos(z4/r4);
    fai4  = atan2(y4,x4);
    if (fai4 < 0)
        fai4 = fai4 + 2*pi;
    end
    ra4 = fai4;
    de4 = pi/2 - sita4;
    % rangeid2
    for i = 1 : gridnum
        if ((ra2>=cat(i).ramin)&(ra2<=cat(i).ramax)&(de2>=cat(i).demin)&(de2<=cat(i).demax))
            rangeid2 = i;
            break;
        end
    end
    % rangeid3
    for i = 1 : gridnum
        if ((ra3>=cat(i).ramin)&(ra3<=cat(i).ramax)&(de3>=cat(i).demin)&(de3<=cat(i).demax))
            rangeid3 = i;
            break;
        end
    end  
    % rangeid4
    for i = 1 : gridnum
        if ((ra4>=cat(i).ramin)&(ra4<=cat(i).ramax)&(de4>=cat(i).demin)&(de4<=cat(i).demax))
            rangeid4 = i;
            break;
        end
    end 
    % catpair
    catpair(igrid1).pair(1) = rangeid1;
    catpair(igrid1).pair(2) = rangeid2;
    catpair(igrid1).pair(3) = rangeid3;  
    catpair(igrid1).pair(4) = rangeid4;
end


