function S = GmKx(quant,kx,f,N,params);
% function S = GmKx(quant,kx,f,N,params);
% Return 1-D horizontal wave-number spectra. 
%
% quant is one of 'U','V','Vel' = U+V,
% 'Uz','Vz','Shear','Disp','Strain'.   
%
% kx are the horizontal wave numbers we want in the spectrum [cpm].  The
% function assumes an isotropic wave field kh=sqrt(kx^2+ky^2).
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
% Optional entries are:
%        .Nphi  - integration resolution in frequency
%        .Nz    - integration resolution in kH
%        .trimlow - 1 implies no wavenumbers corresponding to j<1
%        .trimhigh - 1 implies no wavenumbers above kz=0.1 cpm.  
%        .Nkz  - integration resolution in kz (default=100);

% $Id$
% J. Klymak, April, 2004.  

Nkz=100;
possible={'Nkz'};
for i=1:length(possible);  
  if isfield(params,possible(i));
    eval(sprintf('%s=%f;',possible{i},params.(possible{i})))
  end;
end;  

kz = logspace(-4,1,Nkz);
S = GmKxKz(quant,kx,kz,f,N,params);
S = trapz(kz,S);

