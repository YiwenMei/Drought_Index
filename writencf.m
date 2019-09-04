function writencf(ncn,tag_o,Vb_n,Vb,ndv,Atr_n,Atr,Gatr_n,Gatr)
%% Check for inputs
[~,nm,~]=fileparts(ncn);
fprintf('\nStart processing of %s:%s',nm,Vb_n);
switch nargin
    case 6; error('Not enough number of arguments');

    case 7
%% Variables and variable attributes
      Vb(isnan(Vb))=ndv;
      Vb=Vb';
      switch tag_o
          case 'CW'
            nccreate(ncn,Vb_n,'Dimensions',{'lon',size(Vb,1),'lat',size(Vb,2),'time',size(Vb,3)},...
              'Format','netcdf4_classic','Datatype','double');
            ncwrite(ncn,Vb_n,Vb);
            for a=1:length(Atr_n)
              ncwriteatt(ncn,Vb_n,Atr_n{a},Atr{a});
            end
  
          case 'W'
            ncwrite(ncn,Vb_n,Vb);
            for a=1:length(Atr_n)
              ncwriteatt(ncn,Vb_n,Atr_n{a},Atr{a});
            end

          case 'Loc'
            nccreate(ncn,Vb_n,'Dimensions',{Vb_n,length(Vb)},'Format','netcdf4_classic','Datatype','double');
            ncwrite(ncn,Vb_n,Vb);
            for a=1:length(Atr_n)
               ncwriteatt(ncn,Vb_n,Atr_n{a},Atr{a});
            end

          otherwise; error('Please specify task type');
      end

    case 8; error('Attributes of Geo-information missing');

    case 9
%% Variables and variable attributes
      Vb(isnan(Vb))=ndv;
      Vb=Vb';
      switch tag_o
          case 'CW'
            nccreate(ncn,Vb_n,'Dimensions',{'lon',size(Vb,1),'lat',size(Vb,2),'time',size(Vb,3)},...
              'Format','netcdf4_classic','Datatype','double');
            ncwrite(ncn,Vb_n,Vb);
            for a=1:length(Atr_n)
              ncwriteatt(ncn,Vb_n,Atr_n{a},Atr{a});
            end
  
          case 'W'
            ncwrite(ncn,Vb_n,Vb);
            for a=1:length(Atr_n)
              ncwriteatt(ncn,Vb_n,Atr_n{a},Atr{a});
            end

          otherwise; error('Please specify task type');
      end

%% File attributes
      [~,nm,~]=fileparts(ncn);
      ds=nm(6:end);
      Gatr=[Gatr;datestr(datenum(ds,'yyyymmddHH'),'HH:MM:SS');...
          datestr((24*datenum(ds,'yyyymmddHH')+1)/24,'HH:MM:SS')];

      for a=1:length(Gatr_n)
        ncwriteatt(ncn,'/',Gatr_n{a},Gatr{a});
      end
    
    otherwise; error('Too many number of arguments');
end
fprintf(' ----Done!\n');
end
