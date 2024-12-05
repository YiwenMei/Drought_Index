% Yiwen Mei (meiyw3@mail.sysu.edu.cn), Yuandi Zhou (zhouyd23@mail2.sysu.edu.cn)
% SGP, Sun Yat-sen University
% Last update: 12/4/2024

%% Functionality
% This function is used to calculate the seasonality index. It can be applied
%  to flux or state variables. The user can choose from using either Walsh & Lawler (1981)
%  or Sumner et al. (2001) method for the calculation.

%% Input
% Tbl : Matlab time table stores the monthly flux or state variable time series;
% tflg: Variable type flag which can be 'flux' or 'state' for flux or state variable;
% mflg: Method flag which can be 'S2001' for the Sumner et al. (2001) method
%        or 'WL1981' for the Walsh & Lawler (1981) method.

%% Output
%  SI : Seasonality index of the variable;
% v_Ym: Mean annual accumulated flux or mean state over years.

function [SI,v_Ym]=cal_SI(Tbl,tflg,mflg)
%% Check the inputs
narginchk(3,3);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Tbl',@(x) validateattributes(x,{'timetable'},{'ncols',1},mfilename,'Tbl'));
addRequired(ips,'tflg',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'tflg'));
addRequired(ips,'mflg',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'mflg'));

parse(ips,Tbl,tflg,mflg);
clear ips

%% Fill the incomplete years
Tbl=retime(Tbl,'monthly');
[Y,M,~]=datevec(Tbl.Time);
val_Mm=accumarray(M,Tbl{:,1},[],@(X) mean(X,'omitnan')); % Mean monthly value
k=isnan(Tbl{:,1}); % Index of the missing months
Tbl{k,:}=val_Mm(M(k)); % Fill the missing months by the mean monthly value

yid=Y-min(Y)+1; % Year ID
Nm=accumarray(yid,1); % Number of month of the years
Yid=accumarray(yid,Y,[],@mode); % Unique years
k=Nm<12;
ry=Yid(k); % Incomplete years
kt=ismember(Y,ry);
Tbl{kt,:}=NaN;

%% Calculate the mean and SI
switch tflg
  case 'flux'
    val_Y=accumarray(yid,Tbl{:,1},[],@sum,NaN); % Annual accumulation
    v_Ym=mean(val_Y,'omitnan'); % Mean annual flux

    switch mflg
      case 'S2001'
        Ano=Tbl{:,1}-repelem(val_Y./Nm,Nm,1); % Monthly anormaly
        SI_y=accumarray(yid,abs(Ano),[],@sum)./val_Y;
        SI_y(val_Y==0)=0; % If annual flux/state is 0, SI_y is set to 0
        SI=mean(SI_y,'omitnan'); % Seasonality index

      case 'WL1981'
        val_Mm=accumarray(M,Tbl{:,1},[],@(X) mean(X,'omitnan')); % Mean monthly flux
        SI=sum(abs(val_Mm-v_Ym/12))/v_Ym;
        SI(v_Ym==0)=0;
    end

  case 'state'
    val_Y=accumarray(yid,Tbl{:,1},[],@mean); % mean state
    v_Ym=mean(val_Y,'omitnan'); % mean state over years

    switch mflg
      case 'S2001'
        Ano=Tbl{:,1}-repelem(val_Y,Nm,1); % Monthly anormaly
        SI_y=accumarray(yid,abs(Ano),[],@mean)./val_Y;
        SI_y(val_Y==0)=0; % If annual flux/state is 0, SI_y is set to 0
        SI=mean(SI_y,'omitnan'); % Seasonality index

      case 'WL1981'
        val_Mm=accumarray(M,Tbl{:,1},[],@(X) mean(X,'omitnan')); % Mean monthly state
        SI=mean(abs(val_Mm-v_Ym))/v_Ym;
        SI(v_Ym==0)=0;
    end
end
end
