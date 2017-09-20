clf
fs = 18;
fsbig = 18;
lw = 3;
R = 1;
fslegend = 16;

%pslice = 5;
pslice = 2; 

ha = tight_subplot(2,2,[0.07 0.06],[0.08 0.01],[0.06 0.01]);
%subplot(2,2,3)
axes(ha(3));
good = ~isinf(od.trrs(:,1));
errorbar(od.Dtheta(good,1),od.trrs(good,1),od.trr_errs(good,1),od.trr_errs(good,1),'bo-','LineWidth',lw)
hold on
errorbar(od.Dtheta(good,pslice),od.trrs(good,pslice),od.trr_errs(good,pslice),od.trr_errs(good,pslice),'ks-','LineWidth',lw)
set(gca,'FontSize',fs)
xlabel('D_\psi [unitless]')
ylabel('\tau_r [unitless]')
axis tight

%subplot(2,2,4)
axes(ha(4))
errorbar(od.Dtheta(:,1),od.chisall(:,1),od.chisall_err(:,1),od.chisall_err(:,1),'bo-','LineWidth',lw)
hold on
errorbar(od.Dtheta(:,pslice),od.chisall(:,pslice),od.chisall_err(:,pslice),od.chisall_err(:,pslice),'ks-','LineWidth',lw)
set(gca,'FontSize',fs)
xlabel('D_\psi [unitless]')
ylabel('\chi [unitless]')
axis tight

axes(ha(2));
%subplot(2,2,2)
hold on
hs(1) = errorbar(od.Dtheta(:,1),od.citsall(:,1),od.citseall(:,1),od.citseall(:,1),'bo-','LineWidth',lw);
legs{1} = sprintf('T = %3.3g (simulation)',unique(od.tavgs(:,1)));
%hold on
hs(2) = errorbar(od.Dtheta(:,pslice),od.citsall(:,pslice),od.citseall(:,pslice),od.citseall(:,pslice),'ks-','LineWidth',lw);
legs{2} = sprintf('T = %3.3g (simulation)',unique(od.tavgs(:,pslice)));
hs(3) = plot(od.Dtheta(:,1),od.cosm(:,1),'--','LineWidth',lw,'color',[0 0 0.5]);
legs{3} = sprintf('T = %3.3g (isotropic theory upper bound)',unique(od.tavgs(:,1)));
hold on
hs(4) = plot(od.Dtheta(:,1),od.cosmfull(:,1),'s-.','LineWidth',lw,'color',[0.96 0.63 0.26]);
legs{4} = sprintf('T = %3.3g (full theory)',unique(od.tavgs(:,1)));

hs(5) = plot(od.Dtheta(:,pslice),od.cosm(:,pslice),'--','LineWidth',lw,'color','r');
legs{5} = sprintf('T = %3.3g (isotropic theory upper bound)',unique(od.tavgs(:,pslice)));
hs(6) = plot(od.Dtheta(:,pslice),od.cosmfull(:,pslice),'s-.','LineWidth',lw,'color','g');
legs{6} = sprintf('T = %3.3g (full theory)',unique(od.tavgs(:,pslice)));

%zz=legend(hs,legs);
zz = legend(hs([1 3 4 2 5 6]),legs([1 3 4 2 5 6]));
set(zz,'FontSize',fslegend);
set(gca,'FontSize',fs)
xlabel('D_\psi [unitless]')
ylabel('Chemotactic index CI')
axis tight
ylim([0.4 0.9])


axes(ha(1));
%subplot(2,2,1)
ks = [1 3 6 10]+20;
%oa = 10;
%pos = {[0 0],[oa 0],[0 -oa*1.5],[oa -oa*1.5]};
oa = 8;
%pos = {[0 0],[oa*1.2 0],[oa*2.4 0],[oa*3.6 0]};
pos = {[0 0],[oa*1.2 0],[0 -oa*1.2],[oa*1.2 -oa*1.2]};

for i = 1:length(ks)
    ccv_centered_shot(od.xss{ks(i)}(:,end),od.yss{ks(i)}(:,end),od.thetass{ks(i)}(:,end),od.systematic{ks(i)},od.Sgrad,R,[pos{i} pos{i}]);
    text(pos{i}(1)-oa/3,pos{i}(2)-oa*0.55,sprintf('D_\\psi = %3.3g',od.Dtheta(ks(i))),'FontSize',fsbig);
end
axis equal
axis tight
%xlim([-4 35])
axis off

%ccv_centered_shot(xss,yss,thetass,systematic,Sgrad,R,offset)
