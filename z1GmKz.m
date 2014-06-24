function [S,C] = GmKz(quant,kz,f,N,params);
% function S = GmKz(quant,kz,f,N,params);
% Return vertical wavenumber spectra. 
%
% quant is one of 'U','V','Vel' = U+V,
% 'Uz','Vz','Shear','Disp','Strain'.   
%
% kx are the horizontal wave numbers we want in the spectrum [cpm].
% kz are the vertical wave numbers we want in the spectrum [cpm].
% 
% f is the Coriolis freq [rad/m]
% N is the buoyancy freq [rad/m]
%
% params is a list of parameters.  Must contain:
%  params.s
%        .t 
%        .jstar
%        .jp
%
% For parameters, you can call Gm76Params.m etc...

% $Id$
% J. Klymak, April, 2004.  

b=1300;
N0 = 5.2e-3;
E0 = 6.3e-5;
Nphi = 10000;

possible={'b','N0','E0','Nphi'};
for i=1:length(possible);  
  if isfield(params,possible(i));
    eval(sprintf('%s=%f;',possible{i},params.(possible{i})))
  end;
end;  

% Set the parameters that must be set.  
s=params.s;
t=params.t;
jp = params.jp;
jstar = params.jstar;

% Get normalization constant I
I = s*gamma(t/s)./gamma(1/s)./gamma((t-1)/s);

if jstar<0
  j0=20;jinf=10;
  om = f*sec(phi);
  om0=f;
  ominf = 1.133*2*pi/3600
  omm=0.173*2*pi/3600
  je = j0+0.5*(jinf-j0).*(1-tanh((log10(om/f)-log10(omm/f))./(0.25* ...
    log10(om0/ominf))));
  J=2.1
  jstar = je./J
end;


if size(kz,1)==1
  kz=kz';
end;


C = gmVert(kz,f,jstar,jp,N,b,N0,I,s,t,params);

eps = f/N;
phi = (1:(Nphi))*(acos(eps))/(Nphi+1);
dphi = median(diff(phi));
B = gmFreqPhi(phi,f);
dom =f.*tan(phi).*sec(phi)*dphi;

om = f*sec(phi);

switch upper(quant(1:min(3,end)))
 case upper('dis')
  R = b.^2*N0/N.*(om.^2-f.^2)./om.^2;
 case upper('Vel')
  R = b.^2*N0*N*(om.^2+f.^2)./om.^2;
 case upper('Str')
  R = (2*pi*kz).^2*(b.^2*N0/N.*(om.^2-f.^2)./om.^2);
 case upper('She')
  R = (2*pi*kz).^2*(b.^2*N0*N*(om.^2+f.^2)./om.^2);
 case upper('Ene')
  R = 1+0*kz*om;
 otherwise
  error(sprintf('Do not recognize quant=%s',upper(quant)));
end;

if size(R,1)==1
  R = repmat(R,length(kz),1);
end;

S=trapz((E0*(C*(B.*dom)).*R)');


