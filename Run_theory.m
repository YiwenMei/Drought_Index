% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 8/7/2021

%% Functionality
% This code implement the run theory proposed by Yevjevich (1967).

% Yevjevich, V. M. (1967). Objective approach to definitions and investigations
%  of continental hydrologic droughts, An (Doctoral dissertation, Colorado State
%  University. Libraries).

%% Input
% SDI : Matlab timetable object storing the standardized drought index (SDI);
% thr : the three SDI threshold and the two gap month threshold required by the
%        run theory;
% Nmon: minimum number of month for drought events to be considered.

%% Output
% Rc: all drought events' characteristics stored as a Matlab table object.

%% Additional note
% Requrie run_length.m.

function Rc=Run_theory(SDI,thr,Nmon)
%% Check the inputs
narginchk(3,3);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'SDI',@(x) validateattributes(x,{'timetable'},{'nonempty'},mfilename,'SDI'));
addRequired(ips,'thr',@(x) validateattributes(x,{'double'},{'numel',5},mfilename,'thr'));
addRequired(ips,'Nmon',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Nmon'));

parse(ips,SDI,thr,Nmon);
clear ips

%% Interpolate the SDI to daily
Lo=size(SDI,1);
SDI=retime(SDI,'daily','linear','EndValues',-9999);
Ln=size(SDI,1);
SDI{SDI{:,1}==-9999,1}=NaN;
thr=thr*Ln/Lo;
Nmon=Nmon*Ln/Lo;

%% Identify prelimilary event
k=double(SDI{:,1}<thr(2));

%% Remove insignificant events
[rc,Rc]=Run_Length(k,true,SDI{:,1});
Rc=[reshape(rc',2,length(rc)/2)' Rc'];
Rc(Rc(:,3)>thr(3) & Rc(:,2)<=thr(4),1)=0; % Insignificant event (1-month and higher than -1)

k=Run_Length(reshape(Rc(:,1:2)',size(Rc,1)*2,1),false,[])';

%% Pool together events
[rc,Rc]=Run_Length(k,true,SDI{:,1});
Rc=[reshape(rc',2,length(rc)/2)' Rc'];
Rc(Rc(:,3)<thr(1) & Rc(:,2)<=thr(5),1)=1; % Pool together events separated by 1 month and lower than 1

k=Run_Length(reshape(Rc(:,1:2)',size(Rc,1)*2,1),false,[])';

%% Remove the incomplete events on both ends and short-duration event
[rc,Rc]=Run_Length(k,true,SDI{:,1});
Rc=[reshape(rc',2,length(rc)/2)' Rc'];

if Rc(1,1)>0
  Rc(1,1)=0; % Incomplete event in the beginning
end
if isnan(Rc(1,3)) && Rc(2,1)>0
  Rc(2,1)=0;
end

if Rc(end,1)>0
  Rc(end,1)=0; % Incomplete event in the end
end
if isnan(Rc(end,3)) && Rc(end-1,1)>0
  Rc(end-1,1)=0;
end

Rc(Rc(:,2)<Nmon)=0; % Short-duration event

k=Run_Length(reshape(Rc(:,1:2)',size(Rc,1)*2,1),false,[])';

%% Recode the event mask time series
[rc,Rc]=Run_Length(k,true,SDI{:,1});
Rc=[reshape(rc',2,length(rc)/2)' Rc'];
k1=[-1;diff(Rc(:,1))];
rc=zeros(size(Rc,1),1);
rc(k1==1)=cumsum(Rc(k1==1,1));
Rc(:,1)=rc;

% % Output the event time mask
% TS_mk=Run_Length(reshape(Rc(:,1:2)',size(Rc,1)*2,1),false,[])';
% TS_mk=array2timetable(TS_mk,'RowTimes',SDI.Time,'VariableNames',{'MK_TS'});

%% Event characteristics
Rc=[[1;cumsum(Rc(1:end-1,2))+1] cumsum(Rc(:,2)) nan(size(Rc,1),2) Rc(:,2:3)*Lo/Ln Rc(:,3)./Rc(:,2)...
    nan(size(Rc,1),8)];
Rc(rc==0,:)=[];
for j=1:size(Rc,1)
  T=(Rc(j,1):Rc(j,2))';

  [Sx,Sx_i]=min(SDI{T,1});
%   Rc(j,4)=T(Sx_i)*Lo/Ln; % Time reaches the maximum intensity
  Rc(j,4)=datenum(SDI.Time(T(Sx_i)));
  Rc(j,8)=Sx; % Maximum intensity

  wt=SDI{T,1}/sum(SDI{T,1});
%   Rc(j,3)=sum(wt.*T*Lo/Ln); % Centroid since the beginning of time series
  Rc(j,3)=datenum(SDI.Time(round(sum(wt.*T))));

  Rc(j,9)=(sum(wt.*T)-Rc(j,1)+1)*Lo/Ln; % Centroid since the beginning of event
  Rc(j,10)=(Rc(j,2)-sum(wt.*T)+1)*Lo/Ln; % Centroid since the ending of event
  Rc(j,11)=sum(wt.*T)/mean(T); % Centroid relative to the half-duration
  if any(wt<0)
    Rc(j,12)=sqrt(mean(((T-sum(wt.*T))*Lo/Ln).^2)); % Spreadness
    Rc(j,13)=sqrt(mean((T-sum(wt.*T)).^2))/std(T,1); % Spreadness relative to a time uniform pattern
  else
    Rc(j,12)=std(T*Lo/Ln,wt);
    Rc(j,13)=std(T,wt)/std(T,1);
  end

  Rc(j,14)=(Sx-thr(2))./((T(Sx_i)-Rc(j,1)+1)*Lo/Ln); % Development rate
  Rc(j,15)=(thr(2)-Sx)./((Rc(j,2)-T(Sx_i)+1)*Lo/Ln); % Recovery rate
end

% Rc(:,1:2)=Rc(:,1:2)*Lo/Ln; % Change the time coordinate for begin and end time
Rc(:,1:2)=datenum(SDI.Time(Rc(:,1:2)));

Rc=[Rc Rc(:,6)./Rc(:,12) (Rc(:,6)./Rc(:,12)-thr(2))./Rc(:,9)... % Peakedness, Development gradient
    (thr(2)-Rc(:,6)./Rc(:,12))./Rc(:,10)]; % Recovery gradient

% Output the event characteristics
Rc=array2table(Rc,'VariableNames',{'Tb','Te','Tc','Tx','D','S','Im','Ix','Ct_Tb','Ct_Te','rCt',...
    'Sd','rSd','rst','rcv','Pk','DGrd','RGrd'});
Rc=addvars(Rc,cellfun(@(X) sprintf('Evt_%i',X),num2cell((1:size(Rc,1))'),'UniformOutput',false),...
    'NewVariableNames',{'ETID'},'Before',1);
Rc.Properties.VariableUnits={'-','ith','ith','ith','ith','month','-','-','-','ith_Tb','ith_Te','-',...
    'month','-','1/month','1/month','1/month','1/month^2','1/month^2'};
end
