% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 8/9/2021

%% Functionality
% This code is used to conduct the Mann-Kendall (MK) test for a time series and
%  to calculate the Theil-Sen's (TS) slope of the time series.

%% Input
% TI: time indices of the time series;
% TV: values of the time series;

% alp : significant level for conducting the MK test (default is 0.05);
% Tflg: flag indicating the type of alternative hypothesis to evaluate (i.e.
%        left, right, or both);

%% Output
% H: whether the null hypothesis is rejected (1) or not (0) (null hypothesis - trend
%     absence in the time series).

function [Slp,H,Z,pv]=MK_ST(TI,TV,varargin)
%% Check the inputs
narginchk(2,4);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'TI',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'TI'));
addRequired(ips,'TV',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'TV'));

addOptional(ips,'alp',.05,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'alp'));
addOptional(ips,'Tflg','left',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'Tflg'));

parse(ips,TI,TV,varargin{:});
alp=ips.Results.alp;
Tflg=ips.Results.Tflg;
clear ips varargin

%% Mann-Kendall test
N=length(TV); 
Cid=combnk(1:N,2);

dV=diff(TV(Cid),1,2);
uv=unique(TV(Cid(dV==0))); % Values for tied groups
if ~isempty(uv)
  tp=arrayfun(@(X) length(find(TV==X)),uv); % Count of tied group
  adj=sum(tp.*(tp-1).*(2*tp+5)); % Adjust for tied values
else
  adj=0;
end

S=sum(sign(dV)); % MK statistics
StdS=sqrt((N*(N-1)*(2*N+5)-adj)/18);

% Significant test
if S>=0
  Z=((S-1)/StdS)*(S~=0);
else
  Z=(S+1)/StdS;
end
[H,pv]=ztest(Z,0,1,'alpha',alp,'Tail',Tflg);

%% Theil-Sen's slope estimate
slp=diff(TV(Cid),1,2)./diff(TI(Cid),1,2);
Slp=median(slp);
end
