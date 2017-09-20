rerun = input('Enter 1 to rerun simulation code, 0 to use existing data:')

leftfrac = 0.35;
leftpad = 0.08;
bottomfrac=0.5;
bottompad=0.12;
othermarg = 0.05;

if(rerun)
	KDs = [0.01 0.1 1 10 100]
	sigmads = logspace(-3,0,15)
	qs = 1:4;
	Ns = 1+3*qs+3*(qs.^2)
	nr = 1e5
	gx = 0.05
	gy = 0
	Nits = 2e3
	disp('Running iterations... (may be slow, 5-10 min)')
	[gv,gxm,gym,gall] = full_deltac_deltasig(Nits,Ns,gx,gy,nr,KDs,sigmads);
else
	load('fig1_data.mat')
end

%subplot(2,2,2);
axes('position',[leftfrac+leftpad bottomfrac+bottompad (1-leftfrac-leftpad-othermarg) (1-bottomfrac-bottompad-othermarg)])
g = sqrt(gx^2+gy^2);
plot_sigma_g(Ns,nr,KDs,sigmads,qs,gv,g);

%subplot(2,2,3)
axes('position',[leftpad bottompad leftfrac-leftpad bottomfrac-bottompad]);
chi = chiq(qs);
plot(Ns,(1/g^2)*2*(sigmads(12)^2)./chi,'k--','LineWidth',3);
hold on;
plot(Ns,(1/g^2)*gv(:,12,3),'ko','MarkerSize',12,'LineWidth',3);
set(gca,'xscale','log','yscale','log','FontSize',18);
xlabel('Number of cells'); ylabel('(\sigma_g/g)^2')
%ylim([1e-4 1e-1]/g^2);
ylim([0.1 0.1/g^2]);
xlim([5 65]);
set(gca,'XTick',[7 19 37 61])

Nsf = linspace(10,50,1e3);
hold on
plot(Nsf,(Nsf.^(-2))/(g^2),'r--','LineWidth',2)
text(19,0.2,'\sim N^{-2}','FontSize',24,'color','r')


box on

clearvars -except rerun leftfrac leftpad bottomfrac bottompad othermarg

if(rerun)
	KDs = logspace(-3,3,25);
	sigmads = [0 1e-2 0.1 0.25]; 
	qs = 1;
	Ns = 1+3*qs+3*(qs.^2)
	nr = 1e5
	gx = 0.05
	gy = 0
	Nits = 2e3
	disp('Running iterations... (may be slow, 5-10 min)')
	[gv,gxm,gym,gall] = full_deltac_deltasig(Nits,Ns,gx,gy,nr,KDs,sigmads);

else
	load('fig1b_data.mat')
end

g = sqrt(gx^2+gy^2);

cols = {'r','g','b','k',[0.5 0 0.5]};
syms = {'o','s','d','x','^'};
axes('position',[leftfrac+leftpad bottompad (1-leftfrac-leftpad-othermarg) (bottomfrac-bottompad)])
%subplot(2,2,4)
%clf;
hold on;
chi = chiq(qs);
for k = 1:length(sigmads)
    plot(1./KDs,(1/g^2)*(2/chi)*(sigmads(k).^2+(1/nr)*((KDs+1).^2)./KDs),'k--','LineWidth',3);
    plot(1./KDs,(1/g^2)*squeeze(gv(1,k,:)),syms{k},'color',cols{k},'MarkerSize',12,'LineWidth',3)
    text(0.2,min(gv(1,k,:))*0.7*(1/g^2),sprintf('\\sigma_\\Delta = %3.6g',sigmads(k)),'color',cols{k},'FontSize',18);
end    

set(gca,'xscale','log','yscale','log','FontSize',18);
xlabel('$\bar{c}/K_D$','interpreter','latex'); ylabel('Gradient error (\sigma_g/g)^2');
xlim([1e-3 1e3]);
box on
hold off

	
