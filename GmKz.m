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


Nphi = 1000;

possible={'Nphi'};
for i=1:length(possible);  
  if isfield(params,possible(i));
    eval(sprintf('%s=%f;',possible{i},params.(possible{i})))
  end;
end;  

eps = f/N;
phi = (1:(Nphi))*(acos(eps))/(Nphi+1);
om = f*sec(phi);
% keyboard;
S = GmOmKz(quant,om,kz,f,N,params);
S=trapz(om,S');


