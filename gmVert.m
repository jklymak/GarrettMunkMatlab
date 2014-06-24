function [A] = gMVertKz(kz,f,jstar,jp,N,b,N0,I,s,t,params);
% function [A] = gMVertKz(kz,f,jstar,jp,N,b,N0,I,s,t,varargin);
% return the vertical wavenumber part of the spectrum.  
%
% kz is in cpm, and the output is in cpm^{-1}...
% 
% To specify rolled-off spectra you can specify params.Ef.  If
% params.Ef = 1 the vertical spectrum will roll off at 10 m.  If
% Ef>1 it will roll off earlier.

% $Id: gmVert.m,v 1.1 2008/11/25 22:16:31 jklymak Exp jklymak $ 

% Options...
trimlow=0;
trimhigh=0;
Ef=0;

possible={'trimlow','trimhigh','Ef'};
for i=1:length(possible);  
  if isfield(params,possible(i));
    eval(sprintf('%s=%f;',possible{i},params.(possible{i})))
  end;
end;  

delta = jp*N/N0/2/b;
kstar = jstar*N/N0/2/b;
if size(kstar,2)>1
  [kstar,kz]=meshgrid(kstar,kz);
end
% This is the smallest scale wave possible:
kzmax = 0.1; % Only for GM wave field.  Will need to have an E in here soon.  

A = (1./kstar)*I.*(1+((kz-delta)./kstar).^s).^(-t/s);


if Ef>0 % If this is true, then roll off to k_z^-3 above k_z>0.1 cpm...
  A10 =  (1./kstar)*I.*(1+((0.1-delta)./kstar).^s).^(-t/s);
  Aa = A10*((kz/0.1).^(-3));
  A = min(Aa,Ef*A);
end


if trimlow
  %lowest possible k...
  j0=1;
  k0 = j0*N/N0/2/b;
  bad = find(kz<k0);
  A(bad) = 0;
end;
if trimhigh
  bad = find(kz>kzmax);
  A(bad)=0;
end;
%z=z';