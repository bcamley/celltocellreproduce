rerun = input('Enter 1 to rerun simulation code, 0 to use existing data:')

if(rerun)
disp('Re-running, may take 5-10 minutes...')
SNR = logspace(-1,3,250);
tau = logspace(-1,1,250);
ninterp = 1000; 
np = 1e4;
[V,om,SS,TT,c_interp,sigma_interp] = rigid_rotate_all_redefine(SNR,tau,ninterp,np);

 else
   disp('Loading saved data...')
   load('rigid_rotate_data_broader');
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% do color plot

[c,h] = contour(SS,TT,om,[1e-3 1e-3]); % note we have to do contour at 0+ a small amount 
							      clf;
pcolor(SS,TT,om); shading interp;
set(gca,'xscale','log','yscale','log','FontSize',24);
cb = colorbar;
title(cb,'Rotation speed \omega');
%cb.Label.String=
%cb.Location = 'northoutside';
caxis([0 1]);
hold on
plot(c(1,2:end),c(2,2:end),'k','LineWidth',4);
xlabel('SNR with no rotation','FontSize',24);
ylabel('T / \tau_{rot}','FontSize',24)
SNRlow = logspace(-1,0.5,100);
plot(SNRlow,ones(1,100)/sqrt(2),'y--','LineWidth',5)
text(0.15,0.6/sqrt(2),'Low-SNR_0 limit','color','y','FontSize',24)
SNRhigh = logspace(0.5,3,100);
plot(SNRhigh,sqrt(SNRhigh/2-1/4),'r--','LineWidth',5)
text(20,2,'High-SNR_0 limit','color','r','FontSize',24)

%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Do v(omega) plot
%%
%%%%%%%%%

SNR0=1;

tauex = [0.2 1 5];
lw = [4 6 3];
cols = {'r','k','b'};
stys = {'--','-','-.'};
leglabs = {};
figure

hold on;
for j = 1:length(tauex)
om_max = sqrt(2);
oms = linspace(0,om_max,np);
sigmas = sqrt((1/SNR0)*1./(1+(tauex(j)*oms).^2));
c = cos_average_interp(sigmas,sigma_interp,c_interp);
vmags = real(sqrt(1-0.5*(oms.^2)));
vs = c.*vmags;
hs(j) = plot(oms,vs,stys{j},'color',cols{j},'LineWidth',lw(j));
leglabs{j} = sprintf('T/\\tau_{rot} = %1.1g',tauex(j));
end
set(gca,'FontSize',24);
xlabel('Unitless rotation speed $\omega$','FontSize',24,'interpreter','latex')
ylabel('$\bar{v}_x/v_{\textrm{max}}$','FontSize',24,'interpreter','latex')
legend(hs,leglabs);

%%%%%%%%%
%%
%%% Plot omega as f(R)
%%%%%

Rs = linspace(40,150,300);
rhoc = 3.2e-3;
g = 1e-3;
vmax = 1; % micron/min

sigmad = 0.3;
chi = (pi/4)*rhoc*Rs.^4;
    SNR = (g^2)*chi/(sigmad^2)
      Ts = [20 60 5*60];
figure
hold on
      for j = 1:length(Ts)
		OMR = omega_as_func_R(Rs,Ts(j),sigmad,rhoc,g,vmax,sigma_interp,c_interp,np); % this is the unitful version, with Omega in rad/min

  hrot(j) = plot(Rs,60*OMR,stys{j},'color',cols{j},'LineWidth',lw(j));
legrot{j} = sprintf('T = %d min',Ts(j)); 
end
set(gca,'FontSize',24)
xlabel('Cluster radius [microns]','FontSize',24)
ylabel('Optimal rotation \Omega [rad/hr]','FontSize',24)
legend(hrot,legrot);
