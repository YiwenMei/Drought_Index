% Yiwen Mei (yiwen.mei@uconn.edu)
% CIRCA, University of Connecticut
% Last update: 5/12/2022

%% Functionality
% This code is used to download the USGS gauge data record. The code constructs
%  the downloading link and then use wget for downloading.

%% Input
%  sid : cell array for site ID of the gauge;
% sd/ed: start/end date of the record in yyyy-mm-dd format;
%  pc  : parameter code (00010 - Water temperature (degC), 00045 - Precipitation (in/T);
%         00060 - Discharge (cfs), 00065 - Gage height (ft),
%         00095 - Specific conductivity at 25degC (X10^-6S/cm), 00400 - pH,
%         63680 - Turbidity (FNU)) among others;
%  rc  : time resolution code (iv - original resolution of the measurement,
%         dv - daily aggregation value);
%  pth : name of the output directory.

% pflg: parallel flag (false - default, squential; true - parallel);
% onm : cell array for the output file name.

%% Output
% /<pth>/<pc>-<sid>.txt or /<pth>/<onm>.txt: downloaded files.

function USGS_downloader(sid,sd,ed,pc,rc,pth,varargin)
%% Check the inputs
narginchk(6,8);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'sid',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'sid'));
addRequired(ips,'sd',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'sd'));
addRequired(ips,'ed',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ed'));
addRequired(ips,'pc',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'pc'));
addRequired(ips,'rc',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'rc'));
addRequired(ips,'pth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'pth'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));
addOptional(ips,'onm',cell(length(sid),1),@(x) validateattributes(x,{'cell'},{'nonempty'},...
    mfilename,'onm'));

parse(ips,sid,sd,ed,pc,rc,pth,varargin{:});
pflg=ips.Results.pflg;
onm=ips.Results.onm;
clear ips varargin

%% Downloading the record
url='https://nwis.waterservices.usgs.gov/nwis';
fmt_str='format=rdb';
switch pflg
  case true
    parfor s=1:length(sid)
      USGS_downloader_sub(url,rc,fmt_str,sd,ed,pc,sid{s},pth,onm{s});
%       pause(randi([1 4])*.1); % Pause for some .1secs to behave like a real user
    end

  case false
    for s=1:length(sid)
      USGS_downloader_sub(url,rc,fmt_str,sd,ed,pc,sid{s},pth,onm{s});
    end
end
end

function USGS_downloader_sub(url,rc,fmt_str,sd,ed,pc,sid,pth,onm)
%% Form the link with the inputted options
sd_str=sprintf('startDT=%s',sd);
ed_str=sprintf('endDT=%s',ed);
pc_str=sprintf('parameterCd=%s',pc);
sid_str=sprintf('sites=%s',sid);
url_link=sprintf('%s/%s/?%s&%s&%s&%s&%s',url,rc,fmt_str,sd_str,ed_str,pc_str,sid_str);

%% Downloading
if isempty(onm)
  onm=sprintf('%s-%s',pc,sid);
end
fnm=sprintf('%s.txt',onm);
ofn=fullfile(pth,fnm);
if exist(ofn,'file')~=2
  cmstr=sprintf('wget -nv -O %s "%s"',ofn,url_link);
  [~,wstr]=system(cmstr);
%   if contains(wstr,'[81/81]') % Remove empty file
  if ~isempty(strfind(wstr,'[81/81]'))
    delete(ofn);
  else
    fprintf('%s downloaded\n',fnm);
  end
end
end
