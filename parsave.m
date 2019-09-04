function parsave(filename,variables,vbn,arg1,arg2,arg3)
switch nargin
  case {1, 2}; error('Not enough number of arguments');
  case 3
    eval(sprintf('%s=variables;',vbn));
    save(filename,vbn);
  case 4
    eval(sprintf('%s=variables;',vbn));  
    save(filename,vbn,arg1);
  case 5
    eval(sprintf('%s=variables;',vbn));  
    save(filename,vbn,arg1,arg2);
  case 6
    eval(sprintf('%s=variables;',vbn));  
    save(filename,vbn,arg1,arg2,arg3);

  otherwise; error('Too many number of arguments');
end
end
