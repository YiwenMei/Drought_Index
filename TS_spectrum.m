% Yiwen Mei (yiwen.mei@uconn.edu)
% CIRCA, University of Connecticut
% Last update: 3/19/2022

%% Functionality
% This code implemment the calculation of the time scale spectrum, a Spearman
%  Correlation Coefficient (SCC) matrix on the two time scales of two standardized
%  drought indices (SDI). It also calculate the statistics (mean and standard
%  deviation) of the time scale difference between the two SDIs.

%% Input
% SDI1: Matlab timetable object storing the leading SDI;
% SDI2: Matlab timetable object storing the lagging SDI.

% LT_x: Maximum possible lag time of SDI2 to SDI1 (default is 12, indicating
%        that the response of SDI2 to SDI1 should not exceed 12 time steps);
% alp : Significant level for conducting the test for SCC (default is 0.05).

%% Output
% cc: The SCC matrix on the two SDIs' time scales;
% pv: The p-value for if SCCs are statistically larger than 0 at level of alp;
% Ds: Mean and standard deviation of the time scale difference and the mean SCC.

function [cc,pv,Ds]=TS_spectrum(SDI1,SDI2,varargin)
%% Check the inputs
narginchk(2,4);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'SDI1',@(x) validateattributes(x,{'timetable'},{'nonempty'},mfilename,'SDI1'));
addRequired(ips,'SDI2',@(x) validateattributes(x,{'timetable'},{'nonempty'},mfilename,'SDI2'));

addOptional(ips,'LT_x',12,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'LT_x'));
addOptional(ips,'alp',.05,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'alp'));

parse(ips,SDI1,SDI2,varargin{:});
LT_x=ips.Results.LT_x;
alp=ips.Results.alp;
clear ips varargin

%% Construct the spectrum
[cc,pv]=corr(SDI1,SDI2,'Type','Spearman','Rows','complete','Tail','right');
cc(pv>alp)=NaN;

%% Calculate the weighted mean and standard deviation
TS=1:min(length(SDI2),length(SDI1)-LT_x);

[mxc,mxi]=max(cc,[],1,'omitnan'); % Max MD on HD
if any(~isnan(mxc(TS)))
  wt=mxc(TS)/sum(mxc(TS),'omitnan');
  dp=mxi(TS)-(TS);
  mdP=sum(wt.*dp,'omitnan');
  sdP=std(dp,wt,'omitnan');
  mcc=mean(mxc(TS),'omitnan');
else
  mdP=NaN;
  sdP=NaN;
  mcc=NaN;
end
Ds=[mdP sdP mcc];
end
