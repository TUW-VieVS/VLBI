function nam=sounamivs2iers(ivsnam)

% created by Lucia, 3-9-2014
% modifies by Lucia, 15-10-2014

[iers,ivs]=textread('..\CRF\create\supersource\neededFiles\source_translation_IERS_IVS.dat','%s%s','delimiter','  ','headerlines',4);

 ind = strmatch(ivsnam,iers);
 
 if isempty(ind)
    disp(strcat('source ',ivsnam,' not found in source translation catalogue.'));
    nam=ivsnam;
 else 
    nam=ivs(ind(1)); 
 end
 
 if strcmp(ivsnam,'1253-055')
     nam='3C279   ';
 end
 
