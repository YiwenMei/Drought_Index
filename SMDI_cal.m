% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 8/8/2020

%% Functionality
% This code is used to calculate the Soil Moisture Index (SMI) by Sridhar et al. (2007)
%  and the Soil Moisture Drought Severity (SMDS) by Andreadis et al. (2005).

%% Input
% Obj: a spatiotemporal class object (V2DTCls) for the input soil moisture;
% Nm : number of month to calculate the mean;
% DIT: type of drought index (SMI or SMDS);
% prm: parameters for SMI or SMDS calculations;
% pth: path to store the output drought index files;
% fn : a name for the output drought index files (the output files will have
%       the name fn.Nm.DType.yyyymm.tif in pth);
% ors: coordiante system of the output files (it must be the same as the input.

%% Output
% A list of geotiff files with the name fn.Nm.DType.yyyymm.tif in directory pth.

%% Additional note
% Requrie V2DTCls.m, matV2tif.m.

% Parameter prm for SMI
%  - A 3D array storing two 2D maps with the first one stands for field capacity
%     and the second one for the wilting point of the soil;
%  - A 3-element cell vector storing the calculation details for soil moisture
%     probability (SMP); element 1 is the path to output SMP, element 2 is the
%     name of SMP, and element 3 is a logical indicator (1/0) of whether to caculate
%     SMP or not.

function SMDI_cal(Obj,DIT,Nm,prm,pth,fn,ors)
%% Check the inputs
narginchk(7,7);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Obj',@(x) validateattributes(x,{'V2DTCls'},{'nonempty'},mfilename,'Obj'));
addRequired(ips,'DIT',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'DIT'));
addRequired(ips,'Nm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Nm'));
switch DIT
  case 'SMI'
    [~,~,sz,~]=Obj.GridCls;
    addRequired(ips,'prm',@(x) validateattributes(x,{'double'},{'size',[sz 2]},mfilename,'prm'));
  case 'SMDS'
    addRequired(ips,'prm',@(x) validateattributes(x,{'cell'},{'numel',3},mfilename,'prm'));
  otherwise
    error('Please select from available index');
end
addRequired(ips,'pth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'pth'));
addRequired(ips,'fn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'fn'));
addRequired(ips,'ors',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors'));

parse(ips,Obj,DIT,Nm,prm,pth,fn,ors);
clear ips

%% Time line and spatial info
[Y,~,~]=datevec(Obj.TimeCls('begin'));

mk=~isnan(Obj.readCls(1));
obj.xll=Obj.GIf(1,1)-Obj.GIf(3,1)/2;
obj.yll=Obj.GIf(2,2)-Obj.GIf(3,2)/2;
obj.rs=Obj.GIf(3,1);

for mi=1:12
  switch DIT
    case 'SMI'
%% Soil Moisture Index
% Mean SM
      [SM,ei]=SMDI_cal_sub(mi,Nm,Obj,mk);

% Soil parameters
      if mi==1
        FC=prm(:,:,1);
        FC(~mk)=[];
        WP=prm(:,:,2);
        WP(~mk)=[];
        clear prm
      end

% Drought index
      di=(SM-repmat(WP',1,length(ei)))./repmat((FC-WP)',1,length(ei));
      di=5*di-5; % SWDI: di=10*di;
      ifn=sprintf('%s.%02i',fn,Nm);

% Output SMI
      for y=1:length(ei)
        DI=nan(size(mk)); % Reshape to the land mask
        DI(mk)=di(:,y);
        ofn=fullfile(pth,sprintf('%s.%d%02i.tif',ifn,Y(ei(y)),mi));
        matV2tif(ofn,DI,Obj.ndv,obj,ors,pth)
      end

    case 'SMDS'
%% Soil Moisutre Drought Sevirity
      ifn=prm{2};
      if prm{3}
        fprintf('Execute soil moisture probability calculation\n');
% Reshape SM maps
        [SM,ei]=SMDI_cal_sub(mi,1,Obj,mk);

% Soil moisture probability (SMP)
        di=nan(size(SM));
        for y=1:length(ei)
          di(:,y)=sum(SM<=repmat(SM(:,y),1,length(ei)),2);
        end
        di(isnan(SM))=NaN;
        di=di./(repmat(sum(~isnan(di),2),1,length(ei))+1);

% Output SMP
        for y=1:length(ei)
          DI=nan(size(mk)); % Reshape to the land mask
          DI(mk)=di(:,y);
          ofn=fullfile(prm{1},sprintf('%s.%d%02i.tif',ifn,Y(ei(y)),mi));
          matV2tif(ofn,DI,Obj.ndv,obj,ors,pth)
        end
      end
  end
end

%% Output SMDS
if strcmp(DIT,'SMDS')
% Mean SMP
  nfl=dir(fullfile(prm{1},sprintf('%s.*.tif',ifn)));
  nfl=struct2cell(nfl)';
  nfl=cellfun(@(X) fullfile(prm{1},X),nfl(:,1),'UniformOutput',false);
  Obj=V2DTCls(nfl,'SMP',Obj.ndv,1,0,'Bound',Obj.GIf,0,'begin',NaN, {'yyyymm',...
      sprintf('%s.yyyymm.tif',ifn), 'tif',''});
  for mi=1:12
    [di,ei]=SMDI_cal_sub(mi,Nm,Obj,mk);
    di=1-di;

% Output SMDS
    for y=1:length(ei)
      DI=nan(size(mk)); % Reshape to the land mask
      DI(mk)=di(:,y);
      ofn=fullfile(pth,sprintf('%s.%02i.Empirical.%d%02i.tif',fn,Nm,Y(ei(y)),mi));
      matV2tif(ofn,DI,Obj.ndv,obj,ors,pth)
    end
  end
end
end

function [Cvb,ei]=SMDI_cal_sub(mi,Nm,Obj,mk)
TL=Obj.TimeCls('begin');
[y,m,~]=datevec(TL);
TL=[TL y m];

syi=floor((Nm-2)/12)+2;
ysrt=sort(unique(TL(:,2)));
ei=find(TL(:,3)==mi & TL(:,2)>=ysrt(syi));
si=ei-Nm+1;
ei(si<=0)=[];
si(si<=0)=[];

Cvb=[];
for y=1:length(si)
  Obj1=Obj;
  Obj1.Fnm=Obj.Fnm(si(y):ei(y));
  Sm=[];
  for m=1:length(Obj1.Fnm)
    sm=Obj1.readCls(m);
    sm(~mk)=[];
    Sm=[Sm sm'];
  end
  Cvb=[Cvb nanmean(Sm,2)];
end
end
