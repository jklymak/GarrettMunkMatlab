function p = pdif(f,p);
% function p = pdif(f,p);
% do the first diff for a spectrum

  L=length(f);
  if size(f,1)==L;
    f=f';
  end;
  
  if size(p,1)==L;
    p=p';
    flip = 1;
  else
    flip=0;
  end;  
  
  f = repmat(f,size(p,1),1);
  p = p.*f.^2*4*pi^2;
  
  if flip
    p=p';
  end;
  