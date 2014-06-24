function [A,Z,amax,a0] = GMhorizPhi(Phi,f,jstar,jp,N,b,N0,I,s,t,alphax,Nz,params);
% function A = GMhorizPhi(phi,f,jstar,jp,N,b,N0,z);
% Return the frequency part of the GM spectrum.
%

% $Id$ 


% Options...
trimlow=0;
trimhigh=0;
trimbigzmax=0;
Maxz0=100;
possible={'trimlow','trimhigh','trimbigzmax','Maxz0'};

for i=1:length(possible);  
  if isfield(params,possible(i));
    eval(sprintf('%s=%f;',possible{i},params.(possible{i})))
  end;
end;  



% columns are in phi, rows in z.
%lowest possible z.
j0=1;
a0 = j0/2/N0/b.*f*tan(Phi(1,:));
z0 = sqrt((a0/alphax).^2-1);
bad = find(a0<alphax);
z0(bad)=0;
% try a new method...
if ~trimlow
  z0=0*z0;
end;

if ~trimhigh
  zmax=z0+Maxz0;
  amax = sqrt(zmax.^2+1).*alphax;
else
  % This is the smallest scale wave possible:
  kzmax = 0.1; % Only for GM wave field.  Will need to have an E in here soon.  
  % warning('kzmax set to 0.1 cpm');
  amax = kzmax*f*tan(Phi(1,:))./sqrt(N^2-f^2.*sec(Phi(1,:)).^2);
  fprintf(1,'%e %e %e\n',a0(1),amax(1),alphax);
  zmax = sqrt((amax/alphax).^2-1);
  in = find(zmax.^2<0);
  zmax0=zmax;
  zmax(in)=Maxz0;
  if trimbigzmax
    in = find(zmax-z0>Maxz0);
    zmax(in)=Maxz0+z0(in);
  end;
  if 0
    figure(601);
    subplot(1,2,1);
    semilogy(amax)
    subplot(1,2,2);
    semilogy(zmax)
    pause(0.1);
  end;
end;
zmax = repmat(zmax,size(Phi,1),1);
z0=repmat(z0,size(Phi,1),1);
Z =(0:Nz-1)'*ones(1,size(Phi,2))/(Nz-1);
Z = Z.*(zmax-z0)+z0;
% keyboard;

astar = (f/2/N0/b)*jstar.*tan(Phi);
if jp>0
  delta = (f.*jp/2/N0/b)*tan(Phi);
else
  delta=0;
end;

A = (1./astar)*I.*(1+((alphax.*sqrt(Z.^2+1)-delta)./astar).^s).^(-t/s);
%% set everything that is not possible to 0.   

bad = find(imag(A)>10*eps);
A(bad)=0;
if trimhigh
  allbad = find(amax<alphax);
  A(:,allbad)=0;
end;

if ~isreal(Z)
  % keyboard;
end;

%z=z';