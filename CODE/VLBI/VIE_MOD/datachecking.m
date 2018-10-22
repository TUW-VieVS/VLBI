% 
% ************************************************************************
%   Description:
%   returns values with no repetitions
%  
%   Input:										
%      loading               structure array of loading
% 
%   Output:
%      load_checked        updated structure array with station
%                          correction data
%   Coded: 
%   17 May 2016 by A. Girdiuk
%
%   Revision: 
%  23 Jun 2016 by A. Girdiuk : flags are added
%  26 Jan 2017 by D. Landskron: size of second matrix changed from 3 to 5
%
% *************************************************************************

function [load_checked,flag]=datachecking(loading,mjd1,mjd2)

[mjd_corrected,ia,ic]=unique(loading(:,1));
if size(loading,2)==4
    load_corrected = loading(ia,2:4);
end
if size(loading,2)==5
    load_corrected = loading(ia,2:5);
end
load_checked = [mjd_corrected load_corrected];

      index_mjd1 = load_checked(:,1)<=mjd1;           % logical array 1: before session started
ind_last1 = length(load_checked(index_mjd1,1));	% index in loading array is the closest by time to mjd1

      index_mjd2 = load_checked(:,1)>=mjd2;           % logical array 1: before session finished
ind_last2 = length(load_checked(index_mjd2,1));	% index in loading array is the closest by time to mjd2

if ind_last1>=4 && ind_last2>=4
    flag = 1;
else
    flag = 0;
end
