function S = GmOm(quant,om,f,N,params);
% function S = GmOm(quant,om,f,N,params);
% Return frequency-vertical wavenumber spectra. 
%
% quant is one of 'U','V','Vel' = U+V,
% 'Uz','Vz','Shear','Disp','Strain'.   
%
% om are frequencies we want in the spectrum [rad/s].
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
%        .Nkz resolution of kz integration [10000].  
%
% For parameters, you can call Gm76Params.m etc...
  
% $Id: GmOm.m,v 1.1 2008/11/25 22:18:30 jklymak Exp jklymak $
% J. Klymak, April, 2004.  

if nargin<5
  params= Gm76Params;
end;

Nkz = 1000;

possible={'Nkz'};
for i=1:length(possible);  
  if isfield(params,possible(i));
    eval(sprintf('%s=%f;',possible{i},params.(possible{i})))
  end;
end; 

if size(om,2)==1
  om=om';
end;

kz = logspace(-4,1,Nkz);
% kz = unique([kz linspace(1e-3,1e1,Nkz*0.5)]);
S = GmOmKz(quant,om,kz,f,N,params);
S=trapz(kz,S);
