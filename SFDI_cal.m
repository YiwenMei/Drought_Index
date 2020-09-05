% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 8/12/2020

%% Functionality
% This code is used to calculate standardized streamflow index with different
%  cumulative windows (e.g., SSI1, SSI3, SSI12).

%% Input
%  TL : date of monthly streamflow record in Matlab date number;
%  Q  : monthly streamflow record;
% Mlag: number of month considered to calculate the cumulative streamflow;
% pth : path to store the output drought index files;
%  fn : a name for the output drought index files (the output files will have
%        the name fn.Nm.DType.yyyymm.tif in pth).

% lgflg: flag indicating whether to calculate the drought index based on log-transferred
%         index or not (false - default, original magnitude; true - log streamflow).
% pflg : parallel flag (false - default, squential; true - parallel).

%% Output
% SSI: standardized streamflow index calculated based on the multiple accumulation
%       windows specified in Mlag;

% ofn: file name of the .mat file storing the SSIs.

function [SSI,ofn]=SFDI_cal(TL,Q,Mlag,pth,fn,varargin)
%% Check the inputs
narginchk(5,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'TL',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'TL'));
addRequired(ips,'Q',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'Q'));
addRequired(ips,'Mlag',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'Mlag'));
addRequired(ips,'pth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'pth'));
addRequired(ips,'fn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'fn'));

addOptional(ips,'lgflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'lgflg'));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,TL,Q,Mlag,pth,fn,varargin{:});
lgflg=ips.Results.lgflg;
pflg=ips.Results.pflg;
clear ips varargin

%% Calculate the index
SSI=[];
[y,m,~]=datevec(TL);
TL=[TL y m];
for m=1:length(Mlag)
  Ssi=nan(size(Q));
  switch pflg
    case true
      for mi=1:12
        [ssi,ei,vn]=SFDI_cal_sub(TL,mi,Mlag(m),Q,lgflg);
        Ssi(ei)=ssi;
      end

    case false
      for mi=1:12
        [ssi,ei,vn]=SFDI_cal_sub(TL,mi,Mlag(m),Q,lgflg);
        Ssi(ei)=ssi;
      end
  end

%% Store the index in table
  Tssi=array2table([TL(:,1) Ssi], 'VariableNames', {'Dnum', sprintf('Acc_%02i',Mlag(m))});
  if m==1
    SSI=Tssi;
  else
    SSI=outerjoin(SSI,Tssi,'Keys','Dnum','MergeKeys',1);
  end
end
[y,m,~]=datevec(SSI.Dnum);
y=array2table([y m SSI.Dnum],'VariableNames',{'Year','Month','Dnum'});
SSI=outerjoin(y,SSI,'Keys','Dnum','MergeKeys',1);
SSI=removevars(SSI,'Dnum');

%% Output the SSI
ofn=fullfile(pth,sprintf('%s.%s.mat',vn,fn));
save(ofn,'SSI');
end

function [ssi,ei,vn]=SFDI_cal_sub(TL,mi,mlag,Q,lgflg)
%% Cumulative volume
ei=find(TL(:,3)==mi);
si=ei-mlag+1;
ei(si<=0)=[];
si(si<=0)=[];
V=nan(size(si));
for y=1:length(si)
  V(y)=sum(Q(si(y):ei(y)));
end

% Take the log transform or not
vn='SSI';
switch lgflg
  case true
    if ~isempty(find(V==0, 1))
      V=V+min(Q(Q>0))/100;
    end
    vn=sprintf('%s-ln',vn);
    V=log(V);        
end

%% Standardized streamflow index
ssi=(V-nanmean(V))/std(V,1,'omitnan');
end
