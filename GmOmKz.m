function S = GmOmKz(quant,om,kz,f,N,params);
% function S = GmOmKz(quant,om,kz,f,N,params);
% Return frequency-vertical wavenumber spectra. 
%
% quant is one of 'U','V','Vel' = U+V,
% 'Uz','Vz','Shear','Disp','Strain'.   
%
% om are frequencies we want in the spectrum [rad/s].
% kz are the vertical wave numbers we want in the spectrum [cpm].
% 
% f is the Coriolis freq [rad/s]
% N is the buoyancy freq [rad/s]
%
% params is a list of parameters.  Must contain:
%  params.s
%        .t 
%        .jstar
%        .jp
%
% Optional parameters are:
%        .b thermocline depth scale [1300 m].  
%        .N0 thermocline strat. scale [5.2e-3 rad/s]
%        .E0 Energy level [6.3e-5].
%
% See also: Gm75Params, Gm76Params
  
% $Id$
% J. Klymak, April, 2004.  

b=1300;
N0 = 5.2e-3;
E0 = 6.3e-5;

possible={'b','N0','E0'};
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
if size(om,2)==1;
  om=om';
end;


C = gmVert(kz,f,jstar,jp,N,b,N0,I,s,t,params);
B = gmFreq(om,f);
% B = B./trapz(om,B); 
% this should always integrate to 1 or you will have trouble
% making the 2-D spectra integrable.
switch upper(quant(1:min(3,end)))
 case upper('dis')
  R = b.^2*N0/N.*(om.^2-f.^2)./om.^2;
 case upper('Vel')
  R = b.^2*N0*N*(om.^2+f.^2)./om.^2;
 case upper('Str')
  R = (2*pi*kz).^2*(b.^2*N0/N.*(om.^2-f.^2)./om.^2);
 case upper('She')
  R = (2*pi*kz).^2*(b.^2*N0*N*(om.^2+f.^2)./om.^2);
 otherwise
  error(sprintf('Do not recognize quant=%s',upper(quant)));
end;

if size(R,1)==1
  R = repmat(R,length(kz),1);
end;

S=E0*(C*B).*R;

bad = find(om<f | om>N);
S(:,bad)=0;