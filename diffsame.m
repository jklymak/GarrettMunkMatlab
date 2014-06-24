function x = diffsame(x);
% function x = diffsame(x);
 
if size(x,1)==1 & size(x,2)>1
  x = x';
  flip=1;
else
  flip=0;
end;


x = conv2(diff(x),[1 1]'/2,'full');

if flip
  x=x';
end;
