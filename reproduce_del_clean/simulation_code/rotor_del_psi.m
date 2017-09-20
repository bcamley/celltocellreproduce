function [gxs,gys,xss,yss,vxss,vyss,thetass,phis,systematic,pxests,pyests] = rotor_del_psi(N,dt,steps,L,plotevery,dodraw,tau,svar,kappa,lzero,tavg,Dtheta,Sgrad)
	 EVILFLAG=0; % turn EVILFLAG = 1 to disconnect gradient sensing from motion completely
gxs = zeros(1,steps/plotevery);
gys = zeros(1,steps/plotevery);
pxests = zeros(1,steps/plotevery);
pyests = zeros(1,steps/plotevery);
phis = zeros(1,steps/plotevery); % estimated direction of gradient
Req = 1;
Rdraw = Req;

theta = 2*pi*rand;
[xst,yst] = hex_remove_start(N);
xst = xst-mean(xst);
yst = yst-mean(yst);
R = @(th) [cos(th),-sin(th); sin(th),cos(th)];
zs = R(theta)*[xst ; yst];
xs = zs(1,:);
ys = zs(2,:);
N = length(xs);

psis = 2*pi*rand(size(xs));

systematic = svar*randn(size(xs));
Ss = 1 + systematic + Sgrad*(xs-mean(xs));

tt = linspace(0,2*pi,100);

xss = zeros(N,floor(steps/plotevery));
yss = zeros(N,floor(steps/plotevery));

vxss = zeros(N,floor(steps/plotevery));
vyss = zeros(N,floor(steps/plotevery));

thetass = zeros(N,steps/plotevery);

if(dodraw)
    cmap = jet(1024);
    clf;
    hold on;
    
    for i = 1:N
        hs(i) = plot(xs(i)+(Rdraw/2)*cos(tt),ys(i)+(Rdraw/2)*sin(tt),'r');
        gs(i) = plot(xs(i),ys(i),'.r','MarkerSize',36);
    end
    thetas = psis;
    q = quiver(xs,ys,cos(thetas),sin(thetas),0);
end

range = [min(Ss)*0.6 max(Ss)*1.4];

sdt = sqrt(2*Dtheta*dt);


rrest = rand*2*pi;
pxest = 0.01*cos(rrest);
pyest = 0.01*sin(rrest);

neighs = cell(1,N);

lastwarn('');
for s = 1:steps
    
     DTri = DelaunayTri(xs.',ys.');
     %try
     new_neighs = dt_neighbors(DTri,N); % a lot of this is viciously slow and not very well coded, but doesn't matter much for these small systems
     [msglast,msgidlast] = lastwarn;
     if((strcmp(msgidlast,'MATLAB:singularMatrix')|strfind(msglast,'collinear')|strfind(msglast,'singular')))
             %warning('Had problem with setting neighbors - singular matrix - using old neighbors')
             lastwarn('');
             new_neighs = neighs;
     else
         neighs = new_neighs;
     end
     %catch err
     %    getReport(err)
     %    
     %end
     
     FxI = zeros(size(xs));
     FyI = zeros(size(xs));
     for kk = 1:N
        XIJ = xs(kk)-xs(neighs{kk});
        YIJ = ys(kk)-ys(neighs{kk});
        RIJ = sqrt(XIJ.^2+YIJ.^2);
        FxI(kk) = -kappa*sum((RIJ-lzero).*XIJ./RIJ);            
        FyI(kk) = -kappa*sum((RIJ-lzero).*YIJ./RIJ);
     end

    dxs = xs-mean(xs); dys = ys-mean(ys);
    Ss = 1 + systematic + dxs*Sgrad;
    
    Amat = [ sum(dxs.^2) sum(dxs.*dys) ; sum(dxs.*dys) sum(dys.^2)];
    ghats = linsolve(Amat,[sum(Ss.*dxs) ; sum(Ss.*dys)]);
   
    %ghathatx = ghats(1)/sqrt(ghats(1)^2+ghats(2)^2);
    %ghathaty = ghats(2)/sqrt(ghats(1)^2+ghats(2)^2);  % used to average
    %over the unit vector, but decided that was wrong
    
    pxestn = pxest - (dt/tavg)*(pxest-ghats(1));
    pyestn = pyest - (dt/tavg)*(pyest-ghats(2));
    %size(pxest)

    g_est_angle = atan2(pyest,pxest);

    %pestmag = sqrt(pxest.^2 + pyest.^2);
    %pxesthat = pxest./pestmag;
    %pyesthat = pyest./pestmag;
    %
    %thetasn = thetas + (dt/tau)*real(asin(cos(thetas).*pyesthat - sin(thetas).*pxesthat)) + sdt*randn(size(thetas));
    psisn = psis - (dt/tau)*sin(psis) + sdt*randn(size(psis));  %note sdt has Dtheta in it
    thetas = psis + (1-EVILFLAG)*g_est_angle*ones(size(psis));

    vxs = cos(thetas) + FxI;
    vys = sin(thetas) + FyI;  
    
    xs = xs + vxs*dt;
    ys = ys + vys*dt;
    
    
    
    psis = psisn;
    pxest = pxestn;
    pyest = pyestn;
    
    if(mod(s,plotevery)==0)
        if(dodraw)
            if(isdeployed)
                set(gcf,'Visible','Off');
            end
            xso = xs-L*round(xs/L);
            yso = ys-L*round(ys/L);
            for i = 1:N
                set(hs(i),'XData',xso(i)+(Rdraw/2)*cos(tt));  set(hs(i),'YData',yso(i)+(Rdraw/2)*sin(tt));
                
                colortry = color_interp(Ss(i),cmap,range);
                if(~max(isnan(colortry)))
                    set(gs(i),'XData',xso(i),'YData',yso(i),'Color',colortry);
                else
                    set(gs(i),'XData',xso(i),'YData',yso(i),'Color',[0.7 0.7 0.7]);
                end
            end
            set(q,'XData',xso); set(q,'YData',yso); set(q,'UData',cos(thetas)); set(q,'VData',sin(thetas));
            axis([-L/2 L/2 -L/2 L/2]);
            try
                delete(hh);
            end
            hh=triplot(DTri);
            drawnow;
        end
        
        %printf('Snapshot
        
        xss(:,s/plotevery) = xs;
        yss(:,s/plotevery) = ys;
        vxss(:,s/plotevery) = vxs;
        vyss(:,s/plotevery) = vys;
        thetass(:,s/plotevery) = thetas;
        gxs(s/plotevery) = ghats(1);
        gys(s/plotevery) = ghats(2);
	pxests(s/plotevery) = pxest;
        pyests(s/plotevery) = pyest;
        phis(s/plotevery) = atan2(pyest,pxest);
        
    end
    
    if(N>40 & rem(s,1e3)==0)
        fprintf('Step %d, frac done %3.2g \n',s,s/steps);
    end
    
    
    
end
