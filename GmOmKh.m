function S = GmOmKh(quant,om,kx,f,N,params);
% function S = GmKxKz(quant,kx,kz,f,N,params);
% Return 2-D horizontal-vertical wavenumber spectra. 
%
% quant is one of 'U','V','Vel' = U+V,
% 'Uz','Vz','Shear','Disp','Strain'.   
%
% kx are the horizontal wave numbers we want in the spectrum [cpm].
%    This is not the total kh=sqrt(kx^2+ky^2);
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
% Optional entries are:
%        .Nphi  - integration resolution in frequency
%        .Nz    - integration resolution in kH
%        .trimlow - 1 implies no wavenumbers corresponding to j<1
%        .trimhigh - 1 implies no wavenumbers above kz=0.1 cpm.  

% $Id$
% J. Klymak, April, 2004.  

Nphi = 400;
Nz=500;
trimlow=1;
trimhigh=1;
b=1300;
N0 = 5.2e-3;
E0 = 6.3e-5;


possible={'Nphi','Nz','trimlow','trimhigh','b','N0','E0'};
for i=1:length(possible);  
  if isfield(params,possible(i));
    eval(sprintf('%s=%f;',possible{i},params.(possible{i})))
  end;
end;  

if isfield(params,'freqfunc')
  freqfunc=@params.freqfunc
else
  freqfunc = @gmFreq;
end


% Set the parameters that must be set.  
s=params.s;
t=params.t;
jp = params.jp;
jstar = params.jstar;

if isfield(params,'Ef');
  Ef=params.Ef;
else
  Ef=1;
end;


% Get normalization constant I
I = s*gamma(t/s)./gamma(1/s)./gamma((t-1)/s);

%Nphi = 50;
eps = f/N;
phi = (1:(Nphi))*(acos(eps))/(Nphi+1);
dphi = median(diff(phi));


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



if length(jstar)>1
  jstar = repmat(jstar,Nz,1);
end;

[Om,Kx]=meshgrid(om,kx);
delta = jp*2*N*b*sqrt(Om.^2-f.^2);

astar = jstar/2/N/b*sqrt(Om.^2-f.^2);


S = (1./astar)*I./(1 + ((Kx-delta)./astar).^s).^(t/s);
S = S.*gmFreq(Om,f);

  % Data type dependencies - this is for Z^2.
  switch upper(quant(1:min(3,end)))
   case upper('dis')
    R = b.^2*N0/N.*(Om.^2-f.^2)./Om.^2;
   case upper('Vel')
    R = b.^2*N0*N*(Om.^2+f.^2)./Om.^2;
   case upper('Str')
  Kz = Kx*sqrt((N^2-Om.^2)./(Om.^2-f.^2))
    R = (2*pi*Kz).^2.*b.^2*N0/N.*(Om.^2-f.^2)./Om.^2;
   case upper('She')
  Kz = Kx*sqrt((N^2-Om.^2)./(Om.^2-f.^2))
    R = (2*pi*Kz).^2.*b.^2*N0*N.*(Om.^2+f.^2)./Om.^2;
   otherwise
    error(sprintf('Do not recognize quant=%s',upper(quant)));
  end;
  

%% Some more constants,,,
S= S*2/pi*E0;

S = R.*S;
