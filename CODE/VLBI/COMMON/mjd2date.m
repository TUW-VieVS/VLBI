% This function converts mjd to date
%
% [year] = mjd2date(mjd)
%         to
% [year, month, day, hour, minu, sec] = mjd2date(mjd)
%         or (array based)
% [yr(:,1), mo(:,2), d(:,3), hr(:,4), minu(:,5), sec(:,6)]=mjd2date(mjd(:,1));
%
% ----------------------------------------------
%
% input     mjd [k,1] or [1,k] 
%           varargin1: "noRoundSec" ... if set to _1_, this function does 
%                                       not round to integer seconds
%                      e.g. 59.99999, remains 59.99999
%                                       if set to _2_, this function rounds
%                                       to one decimal point: 59,51 -> 59,5
%                           if not set, or set to 0: 59.51 -> 60
%
% output    yr, mo, d, hr, minu, sec; [k,1] or [1,k] 
%           
%
% not completely tested yet,
% copied mainly from mathworks.com/.../271021
%
%
% created
%   April 1st 2010, Matthias Madzak
%
% changed
%   26 April 2010, mm, changed to array based



% 
 function [year, varargout]=mjd2date(mjd, varargin)
% %% function [year, month, day, hour, minu, sec]=mjd2date(mjd)
% 
% Year required
minArgOut=1; 

%% checks
if nargout<minArgOut
    fprintf('\nNo Output Argument defined\nReturning from mjd2date.m')
end

if  nargout>6
    fprintf('\nToo many output arguments (%1.0d) for mjd2date.m\n--> 6 (year, month, day, hour, min, sec) is maximum.\nAdditional output arguments are set to 0.\n\n', nargout)
end

% varargin
optargin = size(varargin,2);

% if there was 1 (or more) optional input argument
if optargin>0
    if varargin{1}==1
        doNotRound=1;
    elseif varargin{1}==2;
        doNotRound=2;
    else
        doNotRound=0;
    end
else
    doNotRound=0;
end

%% start
mjd=double(mjd);

% preallocation
%sec=zeros(size(mjd));
%min=zeros(size(mjd));
%hour=zeros(size(mjd));

% get hours
hour=floor((mjd-floor(mjd))*24);

% get minutes
minu=floor((((mjd-floor(mjd))*24)-hour)*60);

% get seconds
sec=(((((mjd-floor(mjd))*24)-hour)*60)-minu)*60;

% if varargin{1} not set to one -> round
if doNotRound==0
    sec=round(sec);
elseif doNotRound==2
    % round to first decimal
    sec=round(sec*10)/10;
end
    % OLD: round only when second is very close to integer:
    % sec((sec-floor(sec))>0.99)=round(sec((sec-floor(sec))>0.99));

% change secs and min whose sec==60
minu(sec==60)=minu(sec==60)+1;
sec(sec==60)=0;

% mins==60
hour(minu==60)=hour(minu==60)+1;
minu(minu==60)=0;

% calc jd (yet wrong for hour==24)
jd=mjd+2400000.5;

% if hr==24, correct jd and set hour==0
jd(hour==24)=jd(hour==24)+1;
hour(hour==24)=0;

%% 
% integer julian date
jd_int=floor(jd+0.5);

a=jd_int+32044;
b=floor((4*a+3)/146097);
c=a-floor((b*146097)/4);

d=floor((4*c+3)/1461);
e=c-floor((1461*d)/4);
m=floor((5*e+2)/153);

day=e-floor((153*m+2)/5)+1;
month=m+3-12*floor(m/10);
year=b*100+d-4800+floor(m/10);

%% variable outarg

varargout(1)={month};
varargout(2)={day};
varargout(3)={hour};
varargout(4)={minu};
varargout(5)={sec};

% if more than six output arguments given
if nargout>6
    varargout(6:nargout-minArgOut)={0};
end













% %%
% %% 
% % integer julian date
% jd_int=floor(jd+0.5);
% 
% % Fractional part of day
% jdf=jd-jd_int+0.5;
% 
% % Really next calendar day?
% if jdf>=1
%     jdf=jdf-1;
%     jd_int=jd_int+1;
% end
%     
% %hour=jdf*24; % don't know what that should be for...
% l=jd_int+68569;
% n=floor(4*l/146097);
% 
% l=floor(l)-floor((146097*n+3)/4);
% year=floor(4000*(l+1)/1461001);
% 
% l=l-(floor(1461*year/4))+31;
% month=floor(80*l/2447);
% 
% day=l-floor(2447*month/80);
% 
% l=floor(month/11);
% 
% month=floor(month+2-12*l);
% year=floor(100*(n-49)+year+l);
% 
% 
% 


