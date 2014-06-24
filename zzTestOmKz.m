% Test the omega function.

f = sw_f(30);
N = 5.2e-3;

kz = logspace(-4,-1,100);
om = logspace(log10(f)*1.0001,log10(N),10000);
Ssh = GmKz('Vel',kz,f,N,Gm76Params);
Somkz = GmOmKz('Vel',om,kz,f,N,Gm76Params);



jmkfigure(1,2,0.4);
subplot(1,2,1);
loglog(kz,Ssh);hold on;
loglog(kz,trapz(om,Somkz'),'r');
xlabel('\omega [rad s^{-1}]');
ylabel('\phi_{U_z} [s^-2  (rad s^{-1})^{-1}]');
set(gca,'ylim',[1e-6 1e-3]);

subplot(1,2,2);
loglog(kz,Sstr);
xlabel('\omega [rad s^{-1}]');
ylabel('\phi_{\zeta_z} [(rad s^{-1})^{-1}]');
set(gca,'ylim',[1e-2 1e1]);

if dop
  docprint('OmegaSpec',',TestOm.m','png');
end;

%% Compare shapes (these shouldn't matter)...
jmkfigure(2,2,0.4);clf
S76 = GmOm('Shear',om,f,N,Gm76Params);
S75 = GmOm('Shear',om,f,N,Gm75Params);
S91 = GmOm('Shear',om,f,N,Gk91Params);
SIw = GmOm('Shear',om,f,N,IwexParams);
subplot(1,2,1);
loglog(om,S76,om,S75,om,S91,om,SIw);

S76 = GmOm('Vel',om,f,N,Gm76Params);
S75 = GmOm('Vel',om,f,N,Gm75Params);
S91 = GmOm('Vel',om,f,N,Gk91Params);
SIw = GmOm('Vel',om,f,N,IwexParams);
subplot(1,2,2);
loglog(om,S76,om,S75,om,S91,om,SIw);
