% Yiwen Mei (yiwen.mei@uconn.edu)
% CIRCA, University of Connecticut
% Last update: 4/4/2022

%% Functionality
% This code is used to download the SoilGrids data layer. The code constructs
%  the link and then use the internal websave function or external wget program
%  for downloading.

%% Input
% vb_n : variable name code (see the table below)

%         Organic carbon stock   (t/ha)    - ocs;  pH of water     (\)         - phh2o;
%         Organic carbon content (dg/kg)   - soc;  Nitrogen        (cg/kg)     - nitrogen;
%         Organic carbon density (g/dm^3)  - ocd;  Silt content    (g/kg)      - silt;
%         Clay content           (g/kg)    - clay; Sand content    (g/kg)      - sand;
%         Bulk density           (cg/cm^3) - bdod; Coarse fragment (cm^3/dm^3) - cfvo;
%         Cation exchange capacity at ph 7        (mmol(c))/kg - cec;
%         World Reference Base (2006) Soil Groups (\)          - wrb.

% ver_c: version code (e.g., 2.0.1);
% lyr_c: layer code (0-5cm; 5-15cm; 15-30cm; 30-60cm; 60-100m; 100-200m);
% val_t: value type (mean or uncertainty);
%  lon : longitude range (2 elements from West to East);
%  lat : latitude range (2 elements from South to North);
%  ofn : filename for the data to save to.

% tflg: downloading tool flag (false - default, websave function from matlab;
%        true - wget, make sure to install wget before using in this case).

%% Output
% <ofn>.tif: downloaded file.

%% Additional Note
% Need to make sure wget installed on the computer.

function SoilGrid_downloader(vb_n,ver_c,lyr_c,val_t,lon,lat,ofn,varargin)
%% Check the inputs
narginchk(7,8);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'vb_n',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'vb_n'));
addRequired(ips,'ver_c',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ver_c'));
addRequired(ips,'lyr_c',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'lyr_c'));
addRequired(ips,'val_t',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'val_t'));
addRequired(ips,'lon',@(x) validateattributes(x,{'numeric'},{'numel',2},mfilename,'lon'));
addRequired(ips,'lat',@(x) validateattributes(x,{'numeric'},{'numel',2},mfilename,'lat'));
addRequired(ips,'ofn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ofn'));

addOptional(ips,'tflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'tflg'));

parse(ips,vb_n,ver_c,lyr_c,val_t,lon,lat,ofn,varargin{:});
tflg=ips.Results.tflg;
clear ips varargin

%% Downloading the data
% Downloading options
w_str='https://maps.isric.org/mapserv';
m_str=sprintf('map=/map/%s.map',vb_n);
str1='SERVICE=WCS';
v_str=sprintf('VERSION=%s',ver_c);
str2='REQUEST=GetCoverage';
if strcmp(vb_n,'wrb')
  c_str='COVERAGEID=MostProbable';
elseif strcmp(vb_n,'ocs')
  c_str=sprintf('COVERAGEID=%s_0-30cm_%s',vb_n,val_t);
else
  c_str=sprintf('COVERAGEID=%s_%s_%s',vb_n,lyr_c,val_t);
end
str3='FORMAT=image/tiff';
lon_s=sprintf('SUBSET=long(%.4f,%.4f)',lon(1),lon(2));
lat_s=sprintf('SUBSET=lat(%.4f,%.4f)',lat(1),lat(2));
str4='SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326';
str5='OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326';

% Downloading link
url_link=sprintf('%s?%s&%s&%s&%s&%s&%s&%s&%s&%s&%s',...
    w_str,m_str,str1,v_str,str2,c_str,str3,lon_s,lat_s,str4,str5);

% Execute downloading
if tflg
  cmstr=sprintf('wget -nv -O %s "%s"',ofn,url_link);
  system(cmstr);
else
  websave(ofn,url_link);
  fprintf('%s downloaded\n',ofn);
end
end
