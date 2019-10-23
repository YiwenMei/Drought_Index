% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 12/12/2018

%% Functionality
% This function convert geo-referenced 2D matlab variable to geotif.

%% Input
%  tfn : full name of the output geotiff file;
% matV : input matlab variable;
%  xll : x-coordinate of the lower-left corner of the image;
%  yll : y-coordinate of the lower-left corner of the image;
%  rs  : resolution of the image;
%  ndv : no data value of the image;
%  ors : coordinate system of the image;
% wkpth: working directory of the code.

%% Output
% The output image is tfn.

function matV2tif(tfn,matV,xll,yll,rs,ndv,ors,wkpth)
%% Check inputs
narginchk(8,8);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'tfn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'tfn'));
addRequired(ips,'matV',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'matV'));
addRequired(ips,'xll',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'xll'));
addRequired(ips,'yll',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'yll'));
addRequired(ips,'rs',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'rs'));
addRequired(ips,'ndv',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));
addRequired(ips,'ors',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));

parse(ips,tfn,matV,xll,yll,rs,ndv,ors,wkpth);
clear ips varargin

%% Write Matlab variable to .asc
[~,nm,~]=fileparts(tfn);
afn=fullfile(wkpth,[nm '.asc']);
fid=fopen(afn,'w');
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(size(matV,2))],...
    ['nrows ' num2str(size(matV,1))],['xllcorner ' num2str(xll,12)],['yllcorner '...
    num2str(yll,12)],['cellsize ' num2str(rs)],['NODATA_value ' num2str(ndv)]);
if ~isempty(find(isnan(matV), 1))
  matV(isnan(matV))=ndv; % Assign ndv to NaN
end
dlmwrite(afn,matV,'delimiter',' ','-append');
fclose(fid);

%% Convert .asc to geotiff
fun='gdal_translate';
pr1=sprintf('-a_srs %s',ors);
pr2=sprintf('-a_nodata %i',ndv);
system(sprintf('%s %s %s "%s" "%s"',fun,pr1,pr2,afn,tfn));
delete(afn);
end
