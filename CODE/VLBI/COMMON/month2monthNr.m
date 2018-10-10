% This functiont returns the number of a month (e.g. 3) for the month
% string (e.g. 'March').
%
% INPUT
%       month       string containing the (part of the) month. (e.g.
%                   'March'). Can be string array {'March', 'April'}.
%
% OUTPUT
%       monthNr     array (scalar) containing the number of input months.
%                   Has the same size (column vector or row vector) as the
%                   input).


function [monthNr]=month2monthNr(month)

% if input is char -> make cell
if ischar(month)
    month=cellstr(month);
end

% check if input size is ok
if min(size(month))~=1
    fprintf('Size of input array ''months'' is wrong. Must be vector.\nfile: month2monthNr.m\n');
    return;
end

% make column vector
if size(month,2)>size(month,1)
    month=month';
    sizeEdited=1;
end

% all month numbers
monthNr=zeros(size(month));

% english month names / cell array
allMonthNames={'january January JANUARY jänner Jänner JÄNNER';...
    'february February FEBRUARY februar Februar FEBRUAR';...
    'march March MARCH märz März MÄRZ';...
    'april April APRIL';...
    'may May MAY mai Mai MAI';...
    'june June JUNE juni Juni JUNI';...
    'july July JULY juli Juli JULI';...
    'august August AUGUST';...
    'september September SEPTEMBER';...
    'october October OCTOBER oktober Oktober OKTOBER';...
    'november November NOVEMBER';...
    'december December DECEMBER dezember Dezember DEZEMBER'};

% get number of months
nMonths=size(allMonthNames,1);


% for all input months
for k=1:size(month,1)
    
    %monthNr=~cellfun(@isempty, strfind(allMonthNames, 'jun'));
    for m=1:nMonths
        %a=strfind(allMonthNames{m}, month{k});
        tmp=strfind(allMonthNames{m}, month{k});
        
        if ~isempty(tmp)
            monthNr(k)=m;
            break;
        end
    end
    
    % if current monthNr is still empty -> insert manually
    if monthNr(k)==0
        fprintf('Month ''%s'' not found!\n', month{k});
        inputMonthNr=input('Insert month number: ', 's');
        
        % string 2 numerical
        monthNr(k)=str2double(inputMonthNr);
        
        fprintf('\n');
    end
end

% if size was changed at beginning -> make same vector type (eg row vector)
if exist('sizeEdited', 'var')
    monthNr=monthNr';
end

%fprintf('%0.0f \n', monthNr); 



