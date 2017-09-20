clf
fs = 18;
lw = 3;
R = 1;

oo = 10;

ha = tight_subplot(2,2,[0.07 0.06],[0.08 0.01],[0.06 0.01]);
%subplot(2,2,3)
axes(ha(3));
good = ~isinf(od.trrs(:,1));
errorbar(od.kappas(good,1),od.trrs(good,1),od.trr_errs(good,1),od.trr_errs(good,1),'bo-','LineWidth',lw)
hold on
errorbar(od.kappas(good,oo),od.trrs(good,oo),od.trr_errs(good,oo),od.trr_errs(good,oo),'ks-','LineWidth',lw)
kappafine = linspace(2,7,100);
plot(kappafine,10*exp(0.5*kappafine),'r--','LineWidth',3);
text(5,10*exp(0.5*5)/1.1,'$\sim \exp(\kappa/2)$','FontSize',18,'color','r','interpreter','latex')
set(gca,'FontSize',fs,'yscale','log')
xlabel('\kappa [unitless]')
ylabel('\tau_r [unitless]')

%subplot(2,2,4)
axes(ha(4))
errorbar(od.kappas(:,1),od.chisall(:,1),od.chisall_err(:,1),od.chisall_err(:,1),'bo-','LineWidth',lw)
hold on
errorbar(od.kappas(:,oo),od.chisall(:,oo),od.chisall_err(:,oo),od.chisall_err(:,oo),'ks-','LineWidth',lw)
set(gca,'FontSize',fs)
xlabel('\kappa [unitless]')
ylabel('\chi [unitless]')

axes(ha(2));
%subplot(2,2,2)
hold on
hs(1) = errorbar(od.kappas(:,1),od.citsall(:,1),od.citseall(:,1),od.citseall(:,1),'bo-','LineWidth',lw);
legs{1} = sprintf('T = %3.3g',unique(od.tavgs(:,1)));
%hold on
hs(2) = errorbar(od.kappas(:,oo),od.citsall(:,oo),od.citseall(:,oo),od.citseall(:,oo),'ks-','LineWidth',lw);
legs{2} = sprintf('T = %3.3g',unique(od.tavgs(:,oo)));
hs(3) = plot(od.kappas(:,1),od.cosm(:,1),'--','LineWidth',lw,'color',[0 0 0.5]);

legs{3} = sprintf('T = %3.3g (theory upper bound)',unique(od.tavgs(:,1)));
hold on
plot(od.kappas(:,1),od.cosmfull(:,1),'s-.','LineWidth',lw,'color','g');
hs(4) = plot(od.kappas(:,oo),od.cosm(:,oo),'--','LineWidth',lw,'color',[0.5 0.5 0.5]);
legs{4} = sprintf('T = %3.3g (theory upper bound)',unique(od.tavgs(:,oo)));
plot(od.kappas(:,1),od.cosmfull(:,oo),'s-.','LineWidth',lw,'color','y');
legend(hs,legs)
set(gca,'FontSize',fs)
xlabel('\kappa [unitless]')
ylabel('Chemotactic index CI')


axes(ha(1));
%subplot(2,2,1)
ks = [1 3 6 10];
ks = [1 3 6 10]+90;
oa = 8;
pos = {[0 0],[oa*1.2 0],[0 -oa*1.2],[oa*1.2 -oa*1.2]};

for i = 1:length(ks)
    ccv_centered_shot(od.xss{ks(i)}(:,end),od.yss{ks(i)}(:,end),od.thetass{ks(i)}(:,end),od.systematic{ks(i)},od.Sgrad,R,[pos{i} pos{i}]);
    text(pos{i}(1)-oa/3,pos{i}(2)-oa*0.55,sprintf('\\kappa = %3.3g',od.kappas(ks(i))),'FontSize',18);
end
axis equal
axis tight
axis off

%ccv_centered_shot(xss,yss,thetass,systematic,Sgrad,R,offset)
