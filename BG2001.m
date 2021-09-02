% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 8/1/2020

%% Functionality
% This code perform the change point detection method proposed by Bernaola-Galván
%  et al.(2001).

% Bernaola-Galván, P., Ivanov, P. C., Amaral, L. A. N., & Stanley, H. E. (2001).
%  Scale invariance in the nonstationarity of human heart rate. Phys. Rev. Lett.,
%  87(16), 168105.).

%% Input
%  TS : log-transferred streamflow time series;
%  l0 : the minimum length of a segment where the change point detection will
%        still be performed;
% Nmin: minimum potential change point within a segment with which the change
%        point detection will still be performed;

% alp: significant level for the t-test (default is 0.05).

%% Output
% cpts: the detected change point in time coordinate index;
% TSi : the segment ID series.

function [cpts,TSi]=BG2001(TS,l0,Nmin,varargin)
%% Check the inputs
narginchk(3,4);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'TS',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'TS'));
addRequired(ips,'l0',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'l0'));
addRequired(ips,'Nmin',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Nmin'));

addOptional(ips,'alp',.05,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'alp'));

parse(ips,TS,l0,Nmin,varargin{:});
alp=ips.Results.alp;
clear ips varargin

%% Initialize
cpts=cell(0,0);
flg=true;
TSi=(1:length(TS))';
TS=mat2cell(TS,length(TS),1);
TSi=mat2cell(TSi,length(TSi),1);
Nit=1;
l1=(l0-Nmin)/2;

while flg
  fprintf('The %i-th iteration\n',Nit);

  pts=nan(size(TS));
  pflg=true(size(TS));
  lflg=true(size(TS));
  TS1=cell(2^Nit,1);
  TSi1=cell(2^Nit,1);
  for p=1:length(TS)
    ts=TS{p};
    tsi=TSi{p};
    if sum(~isnan(ts))>=l0
      t_sts=nan(size(ts));

%% Calculate the t statistics
      for t=l1+1:length(ts)-l1
% Statistics of the two segments
        u1=nanmean(ts(1:t-1));
        s1=std(ts(1:t-1),'omitnan');
        N1=sum(~isnan(ts(1:t-1)));
        u2=nanmean(ts(t:length(ts)-l1));
        s2=std(ts(t:length(ts)-l1),'omitnan');
        N2=sum(~isnan(ts(t:length(ts)-l1)));

% t-statistics and significant test 
        S_D=sqrt(((N1-1)*s1^2+(N2-1)*s2^2)/(N1+N2-2)*(1/N1+1/N2));
        t_sts(t)=abs(u1-u2)/S_D;
      end
      [tx,id]=max(t_sts);

%% Detected change points
      N=sum(~isnan(ts));
      eta=4.19*log(N)-11.54;
      v=N-2;
      I=betainc(v/(v+tx^2),.4*v,.4);
      P_tx=(1-I)^eta;

      if P_tx>=alp
        pts(p)=tsi(id);

% Update the time series
        TS1{2*p-1}=ts(1:id);
        TS1{2*p}=ts(id+1:end);
        TSi1{2*p-1}=tsi(1:id);
        TSi1{2*p}=tsi(id+1:end);

% Criteria of ending
      else
        TSi1{2*p-1}=tsi;
        pflg(p)=false;
      end
    else
      TSi1{2*p-1}=tsi;
      lflg(p)=false;
    end
  end
  flg=any(all([pflg lflg],2));

%% Update the time series and record the change points
  if flg
    TS=TS1;
    TSi=TSi1;
    cpts=[cpts pts];
  end
  Nit=Nit+1;
end

% Recode the index time series
a=1;
for i=1:length(TSi)
  if ~isempty(TSi{i})
    TSi{i}=a*ones(size(TSi{i}));
    a=a+1;
  end
end
TSi=cell2mat(TSi);
end
