% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/12/2019

%% Functionality
% This function reads 2D variable stores in different types of file.

%% Input
% Fdetail: details of the 2D variable include the a)name of the file stores the
%          variable, b)its no data value, c)its upper limit, d)its lower limt,
%          and e)its name.

%% Output
% v2d: 2D image as matlab variable.

%% Additional note
% The function supports .tif, .tiff, .nc, .nc4, .hdf, .hdf5, .asc, .txt, and
%  .mat format files.

function v2d=read2Dvar(Fdetail,varargin)
%% Check inputs
narginchk(1,2);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Fdetail',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'Fdetail'));
addOptional(ips,'FunName',mfilename,@(x) validateattributes(x,{'struct'},{'nonempty'},mfilename,'Fdetail'));
parse(ips,Fdetail,varargin{:});

if ~isempty(varargin)
  FunName=ips.Results.FunName;
  FN=[];
  for i=1:length(FunName)
    FN=[FunName(i).name ' -> ' FN];
  end
  FN=[FN mfilename];
else
  FN=mfilename;
end  
clear ips FunName varagin

if iscell(Fdetail)
%% Read the file
  [~,nm,fex]=fileparts(Fdetail{1});
  switch fex
    case {'.tif','tiff'} % compatable for .tiff
      v2d=double(imread(Fdetail{1}));
      nm=[nm fex];
    case {'.nc4','nc'} % compatable for .nc
      v2d=double(ncread(Fdetail{1},Fdetail{5}))';
      nm=[nm fex ':' Fdetail{5}];
    case {'.hdf','hdf5'} % compatable for .hdf5
      v2d=double(hdfread(Fdetail{1},Fdetail{5}));
      nm=[nm fex ':' Fdetail{5}];
    case {'.asc','.txt'}
      v2d=double(dlmread(Fdetail{1},Fdetail{5},Fdetail{6},Fdetail{7}));
      nm=[nm fex ':' Fdetail{5}];
    case '.mat'
      v2d=matfile(Fdetail{1});
      eval(sprintf('v2d=v2d.%s;',Fdetail{5}));
      nm=[nm fex ':' Fdetail{5}];
  end
  v2d(v2d==Fdetail{2})=NaN;

else
  v2d=Fdetail;
  return;
end

%% Check the boundary
validateattributes(v2d(~isnan(v2d)),{'double'},{'<=',Fdetail{3},'>=',Fdetail{4}},'',nm);
% fprintf('%s read by %s\n',nm,FN);
end
