% 
% ************************************************************************
%   Description:
%   subfunction to get corrections for the station of interest
%   takes 4 points before and 4 points after
%  
%   Input:										
%      array               structure array of loading
%      stat                station name of interest
%      mjd                 Modified Julian Date of first observation for 
%                          the scan
%      additional          for 4 values per 1 day
% 
%   Output:
%      mat                 structure array with station
%                          correction data
%   Coded: 
%   17 May 2016 by A. Girdiuk
%
%   Revision: 
%	24 May 2016 by A.Girdiuk: check of data presence
%   23 Jun 2016 by A.Girdiuk: bug-fix
% *************************************************************************
function mat=greploaddata(array,stat,mjd1,mjd2)

index = strcmpi({array.ivsname},stat);	%	 index of station name in loading array

index_mjd1 = array(index).mjd(:,1)<=mjd1;           % logical array 1: before session started
ind_last1 = length(array(index).mjd(index_mjd1,1));	% index in loading array is the closest by time to mjd1

index_mjd2 = array(index).mjd(:,1)<=mjd2;           % logical array 1: before session finished
ind_last2 = length(array(index).mjd(index_mjd2,1));	% index in loading array is the closest by time to mjd2

%index_mjd_2 = array(index).mjd(:,1)>=mjd;			% second posibility is necessary in case of using 2 years of loading data
%ind_last_2 = length(array(index).mjd)-length(array(index).mjd(~index_mjd1,1));

ind_1 = ind_last1-4;				% takes necessary number for spline interpolation
ind_2 = ind_last2+4;            	% and additional for loading which provides data as 4 values per day
									% because of mjd is data of 1st scan,
									% therefore the last scan still needs 4 values for spline interolation
                                    % but old sessions could be longer than 24 hours!!!

if ind_1<1
    ind_1 = 1;							% when we uploaded second year it is true
end
if ind_2>length(array(index).mjd)
    ind_2 = length(array(index).mjd);	% when we need second year it happens
end

mat = array(index).mjd(ind_1:ind_2,:);

