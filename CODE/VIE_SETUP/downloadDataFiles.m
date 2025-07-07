% #########################################################################
% #     downloadDataFiles
% #########################################################################
%
% DESCRITPION
% This function downloads EOP, VMF3, and non-tidal atmospheric loading
% NTAL/VIE files, renames them and places them into the correct folders.
%
% AUTHOR 
%   Sigrid Boehm
%
% INPUT
%   handles     handles.edit_downloadData_period.String string with year or year range
%
% OUTPUT
%
% CHANGES
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function downloadDataFiles(handles)

% Download EOP files independent from dataperiod
try
    c04file = '../EOP/C04_20_1962_now.txt';
    c04url = 'https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now';
    websave(c04file, c04url);
    
    jplfile = '../EOP/JPL_EOP2_long.txt';
    jplurl = 'https://eop2-external.jpl.nasa.gov/eop2/latest_eop2.long';
    websave(jplfile, jplurl);
    
    finfile = '../EOP/finals_all_IAU2000.txt';
    finurl = 'https://maia.usno.navy.mil/ser7/finals2000A.daily.extended';
    websave(finfile, finurl);
    fprintf('EOP files successfully written to ../EOP/.\n');
catch
    warning('There was a problem with your EOP download!');
end

% Check input period and convert to double
dataperiod = handles.edit_downloadData_period.String;
dataperiod = strtrim(dataperiod);

try
    startyr = 0; endyr = 0;
    if length(dataperiod)==4 
        startyr = str2double(dataperiod);
        endyr = startyr;
    elseif length(dataperiod)==9
        startyr = str2double(dataperiod(1:4));
        endyr = str2double(dataperiod(6:9));
    end
    if startyr<1980 || startyr>year(datetime) || endyr<startyr || endyr>year(datetime) || endyr==0
        curyear = dummyfunc(startyr);
    end
    for curyear = startyr:endyr
       vmffile = ['y',num2str(curyear),'.vmf3_r'];
       vmfpath = fullfile('../TRP/VMF3',vmffile);
       if curyear>=1980 && curyear<=2007
           vmfurl = ['https://vmf.geo.tuwien.ac.at/trop_products/VLBI/VMF3/VMF3_EI/yearly/',vmffile];
           websave(vmfpath,vmfurl);
       elseif curyear>=2008 && curyear<=year(datetime)
           vmfurl = ['https://vmf.geo.tuwien.ac.at/trop_products/VLBI/VMF3/VMF3_OP/yearly/',vmffile];
           websave(vmfpath,vmfurl);  
       end
       aplfile = ['vie_y',num2str(curyear),'.ntal_r'];
       aplpath = fullfile('../NTSL/NTAL/VIE',aplfile);
       if curyear>=1981 && curyear<=year(datetime)
           aplurl = ['https://vmf.geo.tuwien.ac.at/APL_products/VLBI/yearly/y',num2str(curyear),'.apl_r'];
           websave(aplpath,aplurl);
       end
       fprintf('VMF3 and NTAL files downloaded successfully.\n');
    end
catch
    warning('off','backtrace');
    warning('There was a problem with your VMF3 and NTAL data download, please check the requested period and try again!');
    warning('on','backtrace');
end