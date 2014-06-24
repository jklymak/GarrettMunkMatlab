% Test the omega function.
dop = 0;

f = sw_f(30);
N = 5.2e-3;

kz = logspace(-4,0,100);
params=Gm76Params;
jmkfigure(1,2,0.4);clf
col = {'r','g','b'}
num = 0;
for Ef = [0 3 1.1];
  num = num+1;
  params.Ef = Ef;
  Ssh = GmKz('Shear',kz,f,N,params);
  Sstr = GmKz('Strain',kz,f,N,params);

  subplot(1,2,1);
  loglog(kz,Ssh,'col',col{num});
  xlabel('k_z [cpm]');
  ylabel('\phi_{U_z} [s^-2  (cpm)^{-1}]');
  set(gca,'ylim',[1e-6 1e-3]);hold on;
  
  subplot(1,2,2);
  loglog(kz,Sstr,'col',col{num});
  xlabel('k_z [cpm]');
  ylabel('\phi_{\zeta_z} [(cpm)^{-1}]');
  set(gca,'ylim',[1e-2 1e1]);hold on;
  
end;
legend('No Roll-off','E=3GM','E=1.1GM',4)
if dop
  print('doc/VerticalSpec','-dpng');
end;

