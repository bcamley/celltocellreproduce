function h = plot_sigma_g(Ns,nr,KDs,sigmads,qs,gv,g)
% function h = plot_sigma_g(Ns,nr,KDs,sigmads,qs,gv,g)
if(nargin<7)
	g = 1;
end
syms = {'o','s','d','x','+'};
cols = {'r','g','b','k'};

chi = chiq(qs);

%clf;
hold on;
% plot theory lines
plot(sigmads,(1/g^2)*2*(sigmads.^2),'r--','LineWidth',3);
text(sigmads(5),(1/g^2)*2*(sigmads(5).^2)/3,'Only cell-cell variation','interpreter','latex','color','r','FontSize',16)

for j = 1:length(KDs)
plot(sigmads,(1/g^2)*2*(sigmads.^2+(1/nr)*((KDs(j)+1)^2)/KDs(j)),'k--','LineWidth',4);
end



for j = 1:length(KDs)
    for k = 1:length(Ns)

    h(j)=plot(sigmads,gv(k,:,j)*chi(k)/g^2,syms{j},'color',cols{k},'MarkerSize',12);

    end
end

for j = 1:length(KDs)

if(j<=ceil(length(KDs)/2))
	text(sigmads(2),1.5*(1/g^2)*2*(sigmads(2).^2+(1/nr)*((KDs(j)+1)^2)/KDs(j)),sprintf('$\\bar{c}/K_D = %2.2g$',KDs(j)),'interpreter','latex','FontSize',16,'color','k')
	end
end

set(gca,'xscale','log','yscale','log','FontSize',18);
xlabel('Cell-to-cell variation \sigma_\Delta','FontSize',18);
ylabel('Rescaled error \chi(\sigma_g/g)^2','FontSize',18)
axis tight
hold off
box on