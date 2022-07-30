function [D,S,I,F]=Drought_Char_LT(SDI,Thr)
for i=1:size(SDI,2)
  sdi=SDI(:,i);
  for t=1:length(Thr)
    thr=Thr(t);
    k=sdi<thr;
    D=sum(k);
    S=sum(sdi(k));
    I=S/D;
    F=D/sum(~isnan(sdi));
  end
end
end
