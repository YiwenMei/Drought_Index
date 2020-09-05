% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 8/1/2020

%% Functionality
% This code is used to calculate water deficit (monthly P-PET). A user supplied
%  PET or PERT calculated by the Thronwite 1948 method is supported.

%% Input
% P_obj: a spatiotemporal class object (V2DTCls) for the input precipitation;
% Mtype: a flag to indicate the type of PET data (T1948 - Thronwite 1948 method,
%         UPET - user-specified PET);
%  obj : a V2DTCls object for the input variable (e.g., monthly air temperature
%         if the T1948 is used, user-specified PET);
% opth : path to store the output drought index files;
%  fn  : a name for the output drought index files (the output files will have
%         the name fn.Nm.DType.yyyymm.tif in pth);

% pflg: parallel flag (false - default, squential; true - parallel).

%% Output
% The monthly water deficit written as geotiff files with the name fn.yyyymm.mat
%  in opth.

%% Additional note
% Requrie date2doy.m and V2DTCls.m.

function WDfc_est(P_obj,Mtype,obj,opth,fn,varargin)
%% Check the inputs
narginchk(5,6);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'P_obj',@(x) validateattributes(x,{'V2DTCls'},{'nonempty'},mfilename,'P_obj'));
addRequired(ips,'Mtype',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'Mtype'));
addRequired(ips,'obj',@(x) validateattributes(x,{'V2DTCls'},{'nonempty'},mfilename,'obj'));
addRequired(ips,'opth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'opth'));
addRequired(ips,'fn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'fn'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,P_obj,Mtype,obj,opth,fn,varargin{:});
pflg=ips.Results.pflg;
clear ips varargin

%% Calculate WDfc
[Y,M,~]=datevec(obj.TimeCls('begin'));
yl=unique(Y);

switch Mtype
%% Thronwite 1984
  case 'T1948'
% Climatological mean temperature
    fprintf('Execute climatological mean air temperature calculation\n');
    rpth=fileparts(obj.Fnm{1});
    for m=1:12
      K=find(M==m);
      N=[];
      T=[];
      for mi=1:length(K)
        T=nansum(cat(3,T,obj.readCls(K(mi))),3);
        n=~isnan(obj.readCls(K(mi)));
        N=sum(cat(3,N,n),3);
      end
      Tclm=T./N;
      ofn=fullfile(rpth,sprintf('Tclm.%02i.mat',m));
      save(ofn,'Tclm');
      fprintf('%s Done\n',ofn);
    end

% Heat index
    fprintf('Execute heat index calculation\n');
    for yi=1:length(yl)
      Hidx=[];
      N=[];
      for m=1:12
        K=find(Y==yl(yi) & M==m);
        if ~isempty(K)
          T=obj.readCls(K);
        else
          ofn=fullfile(rpth,sprintf('Tclm.%02i.mat',m));
          T=matfile(ofn,'T');
          T=T.T;
        end

        I=(T/5).^1.5;
        Hidx=nansum(cat(3,Hidx,I),3);
        n=~isnan(I);
        N=sum(cat(3,N,n),3);
      end
      Hidx=12*Hidx./N;
      ofn=fullfile(rpth,sprintf('Hidx.%d.mat',yl(yi)));
      save(ofn,'Hidx');
      fprintf('%s Done\n',ofn);
    end
    clear T n N Hidx K

    nfl=dir(fullfile(rpth,'Hidx.*.mat'));
    nfl=struct2cell(nfl)';
    nfl=cellfun(@(X) fullfile(rpth,X),nfl(:,1),'UniformOutput',false);
    prm=V2DTCls(nfl,'Hidx',obj.ndv,500,0,'Bound',obj.GIf,0,'begin',NaN,{'yyyy','Hidx.yyyy.mat','mat',...
        ''},'Hidx');

  case 'UPET'
%% User specified PET
    prm=[];
end

fprintf('Execute water deficit calculation\n');
WDfc_est_sub(P_obj,Y,M,Mtype,obj,prm,opth,fn,pflg);
end

function WDfc_est_sub(P_obj,Y,M,Mtype,obj,prm,opth,fn,pflg)
%% Calculate the deficit
switch pflg
  case true
    parfor i=1:length(Y)
      WDfc_est_sub_sub(Mtype,obj,i,Y,M,prm,P_obj,opth,fn);
    end

  case false
    for i=1:length(Y)
      WDfc_est_sub_sub(Mtype,obj,i,Y,M,prm,P_obj,opth,fn);
    end
end
end

function WDfc_est_sub_sub(Mtype,obj,i,Y,M,prm,P_obj,opth,fn)
%% Estimate PET
switch Mtype
  case 'T1948'
    [~,lat,~,~]=obj.GridCls();
    sd=sprintf('%d%02i01',Y(i),M(i));
    nd=eomday(Y(i),M(i));
    ed=sprintf('%d%02i%02i',Y(i),M(i),nd);
    mJD=date2doy((datenum(sd,'yyyymmdd')+datenum(ed,'yyyymmdd'))/2);
    delta=.4093*sin(2*pi()*mJD/365-1.405);
    Nh=24*acos(-tand(lat)*tan(delta))/pi();
    Hidx=prm.readCls(unique(Y)==Y(i));
    alpha=6.75e-7*Hidx.^3-7.71e-5*Hidx.^2+1.79e-2*Hidx+.492;
    T=obj.readCls(i);
    PET=16*(Nh/12).*(nd/30).*(10*T./Hidx).^alpha;

  case 'UPET'
    PET=obj.readCls(i);

  otherwise
    error('Please select from available methods');  
end

%% Water Deficit
WDfc=P_obj.readCls(i)-PET;
WDfc(isnan(WDfc))=obj.ndv;

%% Save the results
ofn=fullfile(opth,sprintf('%s.%i%02i.mat',fn,Y(i),M(i)));
save(ofn,'WDfc');
fprintf('%s Done\n',ofn);
end
