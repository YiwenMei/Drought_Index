% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 12/4/2024

%% Functionality
% This code is used to calculate the Theil-Sen's (TS) slope of a time-depedent 
%  feature and to conduct the Mann-Kendall (MK) test for the TS slope. The unit
%  of the TS slope is in the feature's unit per year.

%% Input
% TI : time indices of the time series (must have at least 3 data points);
% TV : values of the time series (must have at least 3 data points);

% alp: significant level for conducting the MK test (default is 0.05).

%% Output
% Slp: the Theil-Sen's slope;
% Tau: the Kendall Tau correlation coefficient
%  H : whether the null hypothesis is rejected (1) or not (0) (null hypothesis - trend
%       absence in the time series);
% pv : the probability of the Z-Test result.

function [Slp,Tau,H,pv]=TS_MK(TI,TV,varargin)
%% Check the inputs
narginchk(2,3);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'TI',@(x) validateattributes(x,{'datetime'},{'vector'},mfilename,'TI'));
addRequired(ips,'TV',@(x) validateattributes(x,{'double'},{'vector'},mfilename,'TV'));

addOptional(ips,'alp',.05,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'alp'));

parse(ips,TI,TV,varargin{:});
alp=ips.Results.alp;
clear ips varargin

%% Theil-Sen's slope estimate
N=length(TV); 
Cid=nchoosek(1:N,2); % All possible pairs of any two data points

dV=diff(TV(Cid),1,2);
dT_Y=years(diff(TI(Cid),1,2)); % Duration in years
slp=dV./dT_Y; % Unit/years
Slp=median(slp);

%% Mann-Kendall test (equivalent to the Kendall Tau significance test)
Tnum=cumsum([0;days(duration(diff(TI)))]); % Time number
[Tau,pv]=corr(Tnum,TV,'Type','Kendall'); % Double-tail
H=pv<alp;
end
