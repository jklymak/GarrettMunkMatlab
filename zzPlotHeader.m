params= Gm76Params;
kx = logspace(-5,0,200); 
kz = logspace(-4,1,200);
f = sw_f(45);
N=5.2e-3;
om = logspace(log10(f+0.01*f),log10(N),200);
Ef = [0.3 1 3 10];
set(gcf,'defaultlinelinewidth',1.5);

for i=1:length(Ef);
  params.Ef=Ef(i);
  Sz = GmKz('Vel',kz,f,N,params);
  Sx = GmKx('Vel',kx,f,N,params);
  So = GmOm('Vel',om,f,N,params);
  subplot(2,2,1);
  loglog(kz,pdif(kz,Sz));hold on;
  set(gca,'xlim',10.^[-4 1],'xtick',10.^[-5:2]);
  ylabel('S_{dU/dz}(k_z) [s^{-2} cpm^{-1}]');
  xlabel('k_z [cpm]');

  subplot(2,2,2);
  loglog(om,So);hold on;
  set(gca,'xlim',[f N],'xtick',10.^[-5:2]);
  ylabel('S_{U}(\omega) [m^2s^{-2} (rad/s)^{-1}]');
  xlabel('\omega [rad/s]');
  
  subplot(2,1,2);
  loglog(kx,pdif(kx,Sx));hold on;
  set(gca,'xlim',10.^[-5 0],'xtick',10.^[-5:2]);
  ylabel('S_{dU/dx}(k_x) [s^{-2} cpm^{-1}]');
  xlabel('k_x [cpm]');
  
  drawnow
end;