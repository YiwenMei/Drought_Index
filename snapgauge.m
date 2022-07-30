% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 8/20/2021

%% Functionality
% This code performs pour point snapping of stream gauge records based on matching
%  prescribed drainage area with those derived from flow accumulation raster
%  of the domain.

%% Input
% Stn: a Matlab table with three variable under the name of X, Y, and Area representing
%       the horizontal coordinate (m), vertical coordinate (m), and drainage
%       area (km^2);
% fac : the flow accumulation raster (FAC);
% X/Y : the horizontal/vertical coordinate of FAC;
% Nit: maximum number of iteration (every iteration, the maximum possible snapping
%       distance for the gauges is Nit*resolution of FAC;
% Er : tolerance of drainage area (fraction).

%% Output
% stn: the updated station records with the field X, Y, Area, Area_o, and Aflg
%       for the horizontal, vertical coordinate, new, old drainage area, and
%       the area flag indicating good or bad of the matching.

function stn=snapgauge(Stn,fac,X,Y,Nit,Er)
%% Check the inputs
narginchk(6,6);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Stn',@(x) validateattributes(x,{'table'},{'nonempty'},mfilename,'Stn'));
addRequired(ips,'fac',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'fac'));
addRequired(ips,'X',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'X'));
addRequired(ips,'Y',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'Y'));
addRequired(ips,'Nit',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Nit'));
addRequired(ips,'Er',@(x) validateattributes(x,{'double'},{'vector'},mfilename,'Er'));
if length(Er)==1
  Er=Er*ones(Nit,1);
elseif length(Er)~=Nit
  error('Length of Er must be 1 or %i',Nit);
end

parse(ips,Stn,fac,X,Y,Nit,Er);
clear ips

%% Snap the gauges
r_m=mean(abs(diff(X,1,2)),'all');
X=reshape(X,numel(X),1);
Y=reshape(Y,numel(Y),1);

k=true(size(Stn,1),1);
Xn=nan(size(Stn,1),1);
Yn=nan(size(Stn,1),1);
Asn=nan(size(Stn,1),1);
Aflg=cell(size(Stn,1),1);

for i=1:Nit+1
%% Parameters and station inputs
  if i==Nit+1
    Rx=Nit*r_m;
    ErA=1;
  else
    Rx=i*r_m;
    ErA=Er(i);
  end
  fprintf('Iter #%i - ',i);
  fprintf('Search Radius : %im - ',round(Rx));
  fprintf('Tolerence : %.2f%% - ',ErA*100);

  N=round(pi()*(Rx/r_m)^2);
  stn=Stn(k,:);
  Xcor=stn.X;
  Ycor=stn.Y;

%% Perform the kernel calculation
  id=knnsearch([X Y],[Xcor Ycor],'K',N);
  Area=1000^2*stn.Area; % Convert km^2 to m^2
  [dA,j]=min(abs(fac(id)*r_m^2-Area)./Area,[],2,'omitnan');
  I=sub2ind(size(id),(1:length(j))',j);
  I=id(I);
  xn=X(I);
  yn=Y(I);
  An=fac(I)*(r_m/1000)^2; % in km^2
  if i==Nit+1
    aflg=repmat({'bad'},size(stn,1),1);
  else
    aflg=repmat({'good'},size(stn,1),1);
  end

%% Remove overlapped station
  ko=true(size(An));
  [~,m,~]=unique(I,'stable');
  J=I(~logical(accumarray(m,1)));
  for j=1:length(J)
    I_stn=find(I==J(j));
    [~,l]=max(dA(I==J(j)));
    ko(I_stn(l))=false;
  end

%% Output stations that satisfied the criteria
  k1=dA<=ErA & ko; % Station with error of drainage area less than the tolerence
  xn(~k1)=NaN;
  yn(~k1)=NaN;
  An(~k1)=NaN;
  aflg(~k1)={'bad'};
  
  Xn(k)=xn;
  Yn(k)=yn;
  Asn(k)=An;
  Aflg(k)=aflg;

% Update the available pixels
  X(I(k1))=NaN;
  Y(I(k1))=NaN;

% All station satisfied or not
  k=isnan(Asn);
  if sum(k)==0 && i~=Nit+1
    fprintf('%d out of %d stations passed\n',sum(k1),size(stn,1));
    fprintf('All stations passed\n');
    break;
  end
  if i<Nit
    fprintf('%d out of %d stations passed\n',sum(k1),size(stn,1));
  elseif i==Nit
    fprintf('%d out of %d stations passed\n',sum(k1),size(stn,1));
    fprintf('%d out of %d (%.2f%%) stations passed\n',sum(~k),size(Stn,1),100*sum(~k)/size(Stn,1));
  else
    fprintf('Done\n');
  end
end

% Construct the new station records table
stn=array2table([Xn Yn Asn Stn.Area],'VariableNames',{'X','Y','Area','Area_o'});
stn=addvars(stn,(Asn-Stn.Area)./max([Asn Stn.Area],[],2),'NewVariableNames',{'Err_A'});
stn=addvars(stn,Aflg,'NewVariableNames',{'Aflg'});
stn.Properties.VariableUnits={'m','m','km^2','km^2','-','-'};
end
