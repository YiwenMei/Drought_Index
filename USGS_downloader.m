% Yiwen Mei (yiwen.mei@uconn.edu)
% CIRCA, University of Connecticut
% Last update: 2/27/2022

%% Functionality
% This code is used to download the USGS gauge data record. The code constructs
%  the downloading link and then use the websave function for downloading.

%% Input
% sid  : cell array for site ID of the gauge;
% sd/ed: cell array for the start/end date of the record in yyyy-mm-dd format;
%  pc  : cell array for the parameter code (00010 - Water temperature, 00060 - Discharge,
%         00065 - Gage height, 00095 - Specific conductivity at 25degC, 00300 - Dissolved oxygen,
%         00301 - Dissolved oxygen in percent saturation, 00400 - pH);
% pth  : name of the output directory.

% onm : cell array for the output file name;
% pflg: parallel flag (false - default, squential; true - parallel).

%% Output
% /<pth>/<pc>-<sid>.txt or /<pth>/<onm>.txt: downloaded files.

function USGS_downloader(sid,sd,ed,pc,pth,varargin)
%% Check the inputs
narginchk(5,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'sid',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'sid'));
addRequired(ips,'sd',@(x) validateattributes(x,{'cell'},{'numel',length(sid)},mfilename,'sd'));
addRequired(ips,'ed',@(x) validateattributes(x,{'cell'},{'numel',length(sid)},mfilename,'ed'));
addRequired(ips,'pc',@(x) validateattributes(x,{'cell'},{'numel',length(sid)},mfilename,'pc'));
addRequired(ips,'pth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'pth'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));
addOptional(ips,'onm',cell(length(sid),1),@(x) validateattributes(x,{'cell'},{'nonempty'},...
    mfilename,'onm'));

parse(ips,sid,sd,ed,pc,pth,varargin{:});
pflg=ips.Results.pflg;
onm=ips.Results.onm;
clear ips varargin

%% Downloading the record
url='https://nwis.waterservices.usgs.gov/nwis/iv/';
fmt_str='format=rdb';
switch pflg
  case true
    parfor s=1:length(sid)
      USGS_downloader_sub(url,fmt_str,sd{s},ed{s},pc{s},sid{s},pth,onm{s});
    end

  case false
    for s=1:length(sid)
      USGS_downloader_sub(url,fmt_str,sd{s},ed{s},pc{s},sid{s},pth,onm{s});
    end
end
end

function USGS_downloader_sub(url,fmt_str,sd,ed,pc,sid,pth,onm)
%% Form the link with the inputted options
sd_str=sprintf('startDT=%s',sd);
ed_str=sprintf('endDT=%s',ed);
pc_str=sprintf('parameterCd=%s',pc);
sid_str=sprintf('sites=%s',sid);
url_link=sprintf('%s?%s&%s&%s&%s&%s',url,fmt_str,sd_str,ed_str,pc_str,sid_str);

%% Downloading
if isempty(onm)
  onm=sprintf('%s-%s.txt',pc,sid);
end
ofn=fullfile(pth,onm);
if exist(ofn,'file')~=2
  websave(ofn,url_link);
  fprintf('%s downloaded\n',ofn);
else
  fprintf('%s already exists\n',ofn);
end
end
