function [Output,O2]=Run_Length(Input,iflg,TS)
j=1;
k=1;
i=1;

if iflg
%% Run length
  while i<2*length(Input)
    comp=1;
    cumv=0;
    for j=j:length(Input)
      if j==length(Input)
        cumv=cumv+TS(j);
        break;
      end
      if Input(j)==Input(j+1)
        comp=comp+1;
        cumv=cumv+TS(j);
      else
        cumv=cumv+TS(j);
        break;
      end
    end
    Output(k+1)=comp;
    Output(k)=Input(j);
    O2((k+1)/2)=cumv;

    if j==length(Input) && Input(j-1)==Input(j)
      break;
    end
    i=i+1;
    k=k+2;
    j=j+1;
    if j==length(Input)
      if mod(length(Input),2)==0 
        Output(k+1)=1;
        Output(k)=Input(j);
        O2((k+1)/2)=cumv;
      else
        Output(k+1)=1;    
        Output(k)=Input(j);
        O2((k+1)/2)=cumv;
      end

      break;
    end
  end

else
%% Reverse run length
  while i<=length(Input)
    while j<=Input(i+1)
      Output(k)=Input(i);
      j=j+1;
      k=k+1;
    end
    i=i+2;
    j=1;
  end
  O2=[];
end
end
