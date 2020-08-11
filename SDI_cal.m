% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 7/25/2020

%% Functionality
% This code is used to calculate standardized drought index (e.g., SPI, SPEI).

%% Input
%  Obj : a spatiotemporal class object (V2DTCls) for the input variable (e.g.,
%         monthly precipitation, monthly water deficit);
%  Nm  : number of month to calculate the accumulations;
% DType: distribution assumed to the cummulative variable for calculation of
%         the drought index (see Additional note section for distribution available);
%  pth : path to store the output drought index files;
%  fn  : a name for the output drought index files (the output files will have
%         the name fn.Nm.DType.yyyymm.tif in pth);
%  ors : coordiante system of the output files (it must be the same as the input.

%% Output
% A list of geotiff files with the name fn.Nm.DType.yyyymm.tif in directory pth.

%% Additional note
% Requrie matV2tif.m.

% Available distribution for drought index calculation:
%  'Empirical'    - empirical CDF
%  'Exponential'  - exponential distribution
%  'Gamma'        - gamma distribution
%  'Gaussian'     - normal distribution
%  'Log-normal'   - log normal distribution
%  'Log-logistic' - log logistic distribution
%  'Weibull'      - Weibull distribution

function SDI_cal(Obj,Nm,DType,pth,fn,ors)
%% Check the inputs
narginchk(6,6);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Obj',@(x) validateattributes(x,{'V2DTCls'},{'nonempty'},mfilename,'Obj'));
addRequired(ips,'Nm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Nm'));
addRequired(ips,'DType',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'DType'));
addRequired(ips,'pth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'pth'));
addRequired(ips,'fn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'fn'));
addRequired(ips,'ors',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors'));

parse(ips,Obj,Nm,DType,pth,fn,ors);
clear ips varargin

%% Time line and spatial info
TL=Obj.TimeCls('begin');
[y,m,~]=datevec(TL);
TL=[TL y m];
syi=floor((Nm-2)/12)+2;
ysrt=sort(unique(TL(:,2)));

mk=Obj.readCls(1);
mk=~isnan(mk);
xll=Obj.GIf(1,1)-Obj.GIf(3,1)/2;
yll=Obj.GIf(2,2)-Obj.GIf(3,2)/2;

for mi=1:12
%% Cumulative variable
  ei=find(TL(:,3)==mi & TL(:,2)>=ysrt(syi));
  si=ei-Nm+1;
  ei(si<=0)=[];
  si(si<=0)=[];
  Cvb=[];
  for y=1:length(si)
    Obj1=Obj;
    Obj1.Fnm=Obj.Fnm(si(y):ei(y));
    FX=[];
    for m=1:length(Obj1.Fnm)
      vb=Obj1.readCls(m);
      vb(~mk)=[];
      FX=[FX vb'];
    end
    Cvb=[Cvb sum(FX,2)];
  end
  clear Obj1 vb si cvb

%% Calculate the drought index
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
        prm=expfit(Cvb(i,:));
        FX(i,:)=expcdf(Cvb(i,:),prm);
      end

    case 'Gamma'
      for i=1:size(Cvb,1)
        X=Cvb(i,:);
        prm=gamfit(X(X~=0));
        q=sum(X==0)/length(X);
        FX(i,:)=q+(1-q)*gamcdf(X,prm(1),prm(2));
      end

    case 'Gaussian'
      for i=1:size(Cvb,1)
        [prm1,prm2]=normfit(Cvb(i,:));
        FX(i,:)=normcdf(Cvb(i,:),prm1,prm2);
      end

    case 'Log-normal'
      for i=1:size(Cvb,1)
        X=Cvb(i,:);
        prm=lognfit(X(X~=0));
        FX(i,:)=logncdf(Cvb(i,:),prm(1),prm(2));
      end

    case 'Log-logistic'
      Cvb=Cvb-repmat(min(Cvb,[],2),1,size(Cvb,2))+1;
      for i=1:size(Cvb,1)
        prm=fitdist(Cvb(i,:)','loglogistic');
        FX(i,:)=cdf(prm,Cvb(i,:));
      end

    case 'Weibull'
      for i=1:size(Cvb,1)
        prm=wblfit(Cvb(i,:));
        FX(i,:)=wblcdf(Cvb(i,:),prm(1),prm(2));
      end

    otherwise
      error('Please choose from the available list of distribution');
  end
  FX=norminv(FX); % Drought index
  clear Cvb prm X

%% Output the index
  for y=1:length(ei)
    SDI=nan(size(mk)); % Reshape to the land mask (recycle Cvb)
    SDI(mk)=FX(:,y);
    ofn=fullfile(pth,sprintf('%s.%02i.%s.%d%02i.tif',fn,Nm,DType,TL(ei(y),2),mi));
    matV2tif(ofn,SDI,xll,yll,Obj.GIf(3,1),Obj.ndv,ors,pth);
  end
end
end
