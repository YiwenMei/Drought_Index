% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 1/15/2021

%% Functionality
% This code is used to calculate standardized drought index (e.g., SPI, SPEI).

%% Input
%  Obj : a spatiotemporal class object (V2DTCls) or a Matlab timetable for the
%         input variable (e.g. monthly precipitation, monthly water deficit);
%  Nm  : number of month to calculate the accumulations;
% DType: distribution assumed to the cummulative variable for calculation of
%         the drought index (see Additional note section for distribution available);
%  pth : path to store the output drought index files;
%  fn  : a name for the output drought index files (the output files will have
%         the name fn.Nm.DType.yyyymm.tif in pth);

% pflg: parallel flag (false - default, squential; true - parallel).

%% Output
% If Obj is V2DTCls, the 2D drought index for every time step is outputted as
%  geotiff files with the name fn.Nm.DType.yyyymm.tif in pth;
% If Obj is table, the drought index time series for every station is outputted
%  as Matlab table with the name fn.Nm.DType.Sid.mat in pth.

%% Additional note
% Requrie matV2tif.m and V2DTCls.m.

% Available distribution for drought index calculation:
%  'Empirical'    - empirical CDF
%  'Exponential'  - exponential distribution
%  'Gamma'        - gamma distribution
%  'Gaussian'     - normal distribution
%  'Log-normal'   - log normal distribution
%  'Log-logistic' - log logistic distribution
%  'Weibull'      - Weibull distribution

function SDI_cal(Obj,Nm,DType,pth,fn,varargin)
%% Check the inputs
narginchk(5,6);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Obj',@(x) validateattributes(x,{'V2DTCls','timetable'},{'nonempty'},mfilename,'Obj'));
addRequired(ips,'Nm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Nm'));
addRequired(ips,'DType',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'DType'));
addRequired(ips,'pth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'pth'));
addRequired(ips,'fn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'fn'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,Obj,Nm,DType,pth,fn,varargin{:});
pflg=ips.Results.pflg;
clear ips varargin

%% Time line and spatial info
if isa(Obj,'V2DTCls')
  mk=Obj.readCls(1);
  mk=~isnan(mk);
  xll=Obj.GIf(1,1)-Obj.GIf(3,1)/2;
  yll=Obj.GIf(2,2)-Obj.GIf(3,2)/2;

  TL=Obj.TimeCls('begin');

elseif isa(Obj,'timetable')
  mk=1;
  xll=[];
  yll=[];

  TL=datenum(Obj.Time);
  Obj=Obj(:,1:end);  
end
[y,m,~]=datevec(TL);
TL=[TL y m];

%% Calculate the index
ofl=[];
switch pflg
  case true
    parfor mi=1:12
      ofn=SDI_cal_sub(TL,mi,Nm,Obj,mk,DType,pth,fn,xll,yll);
      ofl=[ofl;ofn];
    end

  case false
    for mi=1:12
      ofn=SDI_cal_sub(TL,mi,Nm,Obj,mk,DType,pth,fn,xll,yll);
      ofl=[ofl;ofn];
    end
end
clear TL y mk xll yll mi Nm

%% Concatenate the time series for stations
if isa(Obj,'timetable')
  SDI=[];
  for m=1:12
    sdi=matfile(ofl(m,:));
    sdi=sdi.SDI;
    SDI=[SDI;sdi];

%     delete(ofl(m,:));
  end
  SDI=sortrows(SDI,'Time','ascend');
  clear sdi Obj

% Output the resulting TS
  stn=SDI.Properties.VariableNames;
  for s=1:length(stn)
    ofn=[ofl(1,1:length(pth)+length(fn)+length(DType)+6) stn{s} ofl(1,end-3:end)];
    sdi=SDI(:,s);
    save(ofn,'sdi');
    [~,ofn,~]=fileparts(ofn);
    fprintf('%d. %s saved\n',s,ofn);
  end
end
end

function ofn=SDI_cal_sub(TL,mi,Nm,Obj,mk,DType,pth,fn,xll,yll)
%% Cumulative variable
fprintf('Execute variable accumulation');
ei=find(TL(:,3)==mi);
si=ei-Nm+1;
ei(si<=0)=[];
si(si<=0)=[];
Cvb=[];
for y=1:length(si)
  if isa(Obj,'V2DTCls')
    Obj1=Obj;
    Obj1.Fnm=Obj.Fnm(si(y):ei(y));
    FX=[];
    for m=1:length(Obj1.Fnm)
      vb=Obj1.readCls(m); % vb is 2D
      vb(~mk)=[]; % flatten vb to 1-by-Nloc
      FX=[FX vb']; % flip vb to Nloc-by-1, FX is Nloc-by-Mmonth
    end

  elseif isa(Obj,'timetable')
    FX=Obj{si(y):ei(y),:}';
  end
  Cvb=[Cvb sum(FX,2)]; % sum over month
end
clear Obj1 vb si cvb
fprintf(' ---- Done\n');

%% Calculate the drought index
fprintf('Execute DI calculation - assume %s distribution',DType);
FX=nan(size(Cvb));
switch DType
  case 'Empirical'
    for y=1:length(ei)
      FX(:,y)=sum(Cvb<=repmat(Cvb(:,y),1,length(ei)),2);
    end
    FX(isnan(Cvb))=NaN;
    FX=FX./(repmat(sum(~isnan(Cvb),2),1,length(ei))+1); % Probability

  case 'Exponential'
    for i=1:size(Cvb,1)
      X=Cvb(i,:);
      k=isnan(X);
      X(k)=[];
      prm=expfit(X);
      FX(i,~k)=expcdf(X,prm);
    end

  case 'Gamma'
    for i=1:size(Cvb,1)
      X=Cvb(i,:);
      k=isnan(X);
      X(k)=[];
      prm=gamfit(X(X~=0));
      q=sum(X==0)/length(X);
      FX(i,~k)=q+(1-q)*gamcdf(X,prm(1),prm(2));
    end

  case 'Gaussian'
    for i=1:size(Cvb,1)
      X=Cvb(i,:);
      k=isnan(X);
      X(k)=[];
      [prm1,prm2]=normfit(X);
      FX(i,~k)=normcdf(X,prm1,prm2);
    end

  case 'Log-normal'
    for i=1:size(Cvb,1)
      X=Cvb(i,:);
      k=isnan(X);
      X(k)=[];
      prm=lognfit(X(X~=0));
      q=sum(X==0)/length(X);
      FX(i,~k)=q+(1-q)*logncdf(X,prm(1),prm(2));
    end

  case 'Log-logistic'
    Cvb=Cvb-repmat(min(Cvb,[],2),1,size(Cvb,2))+1;
    for i=1:size(Cvb,1)
      X=Cvb(i,:);
      k=isnan(X);
      X(k)=[];
      prm=fitdist(X','loglogistic');
      FX(i,~k)=cdf(prm,X);
    end

  case 'Weibull'
    for i=1:size(Cvb,1)
      X=Cvb(i,:);
      k=isnan(X);
      X(k)=[];
      prm=wblfit(X(X~=0));
      q=sum(X==0)/length(X);
      FX(i,~k)=q+(1-q)*wblcdf(X,prm(1),prm(2));
    end

  otherwise
    error('Please choose from the available list of distribution');
end
FX=norminv(FX); % Drought index
clear Cvb prm X
fprintf(' ---- Done\n');

%% Output the index
if isa(Obj,'V2DTCls')
  for y=1:length(ei)
    SDI=nan(size(mk)); % Reshape to the land mask (recycle Cvb)
    SDI(mk)=FX(:,y);
    ofn=fullfile(pth,sprintf('%s.%02i.%s.%d%02i.tif',fn,Nm,DType,TL(ei(y),2),mi));
    matV2tif(ofn,SDI,xll,yll,Obj.GIf(3,1),Obj.ndv,Obj.srs,pth);
%     fprintf('Output %s\n',ofn);
  end
  ofn=[];

elseif isa(Obj,'timetable')
  ofn=fullfile(pth,sprintf('%s.%02i.%s.%02i.mat',fn,Nm,DType,mi));
  SDI=array2timetable(FX','RowTimes',datetime(TL(ei,2),TL(ei,3),1),...
      'VariableNames',Obj.Properties.VariableNames);
  save(ofn,'SDI');
end
end
