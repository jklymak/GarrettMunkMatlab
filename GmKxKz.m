function S = GmKxKz(quant,kx,kz,f,N,params);
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

S = ones(length(kz),length(kx));

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


% d omega
dom = repmat(f.*tan(phi).*sec(phi)*dphi,Nz,1);
Phi = repmat(phi,Nz,1);
if length(jstar)>1
  jstar = repmat(jstar,Nz,1);
end;

for i=1:length(kx);
  
  Z = repmat((0:Nz-1)',1,length(kz))./(Nz-1);
  Kz=repmat(kz,size(Z,1),1);
  zmax = sqrt(Kz.^2./kx(i).^2-1);
  Z = Z.*zmax;
  omsq = kx(i).^2./Kz.^2.*(Z.^2+1).*(N.^2-f.^2) + f.^2;
  om = sqrt(omsq);
  if jstar<0
    j0=20;jinf=10;
    om0=f;
    ominf = 1.133*2*pi/3600;
    omm=0.173*2*pi/3600;
    je = j0+0.5*(jinf-j0).*(1-tanh((log10(om/f)-log10(omm/f))./(0.25* ...
      log10(om0/ominf))));
    J=2.1;
    jstar = je./J;
  end;
  if Ef>1
    params.Ef=Ef;
  end;
  C = gmVert(Kz,f,jstar,jp,N,b,N0,I,s,t,params);
  
  B = freqfunc(om,f);
  % dom/da
  domda = kx(i).*sqrt(Z.^2+1)./om./Kz.^2.*(N.^2-f.^2);

  dz= repmat(diff(Z(1:2,:)),size(Z,1),1);  
  % Tda cancels stuff, so just do that here and save some time...
  Tda = 1./sqrt(Z.^2+1).*dz;  

  % Data type dependencies - this is for Z^2.
  switch upper(quant(1:min(3,end)))
   case upper('dis')
    R = b.^2*N0/N.*(om.^2-f.^2)./om.^2;
   case upper('Vel')
    R = b.^2*N0*N*(om.^2+f.^2)./om.^2;
   case upper('Str')
    R = (2*pi*Kz).^2.*b.^2*N0/N.*(om.^2-f.^2)./om.^2;
   case upper('She')
    R = (2*pi*Kz).^2.*b.^2*N0*N.*(om.^2+f.^2)./om.^2;
   otherwise
    error(sprintf('Do not recognize quant=%s',upper(quant)));
  end;
  
  TT = B.*R.*C.*Tda.*domda;
  S(:,i) = trapz(TT)';
end;

%% Some more constants,,,
S= S*2/pi*E0;

S(abs(imag(S))>0)=0;

