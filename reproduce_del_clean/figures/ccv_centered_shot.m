function q = ccv_centered_shot(xss,yss,thetass,systematic,Sgrad,R,offset)
% agh, note R is diameter, not radius... 
N = size(xss,1);

dxss = xss-ones(N,1)*mean(xss);
dyss = yss-ones(N,1)*mean(yss);
Ss = dxss*Sgrad + (systematic.')*ones(1,size(xss,2));
Smax = max(Ss(:));
Smin = min(Ss(:));

tt = linspace(0,2*pi,100);

%cmap = jet(1024);
cmap = parula(1024);

for s = 1:size(xss,2)
    for k = 1:size(xss,1)
       cols{k} = color_interp(Ss(k,s),cmap,[Smin Smax]);
    end  
   %clf;
   hold on;
   xs = dxss(:,s)+offset(1); ys = dyss(:,s)+offset(2);  % always center
   DTri = DelaunayTri(xs,ys);
   hh=triplot(DTri,'color',[0.8 0.8 0.8]);
   set(hh,'linewidth',2);
   %xs = xs-xm*round(xs/xm);
   %ys = ys-ym*round(ys/ym);
   
   q = quiver(xs,ys,cos(thetass(:,s)),sin(thetass(:,s)),0);
   set(q,'LineWidth',2,'color',[0.7 0.7 0.7]);
   
   for k = 1:size(xss,1)
       patch((R/2)*cos(tt)+xs(k),(R/2)*sin(tt)+ys(k),cols{k},'edgecolor','none','facealpha',0.75);
       
   end
   
   
   
   %axis([-L/2 L/2 -L/2 L/2]);
   % axis equal
   %axis([-Lx/2 Lx/2 -Ly/2 Ly/2]);
   %set(gca,'FontSize',18);
   %box on;

   % drawnow
    
end