function untarVgosDb( curNcFolder )
% untar the vgosDb session folder specified with curNcFolder
%
% Input:
%   curNcFolder
% Output:
%   -

vgosTargz = [curNcFolder(1:end-1),'.tar.gz'];
vgosTgz = [curNcFolder(1:end-1),'.tgz'];
curSlash = sort([strfind(curNcFolder,'/'), strfind(curNcFolder,'\')]);
vgosTgzFolder = curNcFolder(1:curSlash(end-1));

if exist(vgosTgz,'file')
    untar(vgosTgz,vgosTgzFolder);
elseif exist(vgosTargz,'file')
    untar(vgosTargz,vgosTgzFolder);
else
    fprintf('ERROR: %s does not exist!\n',vgosTgz);
end

end

