function [vxs,vys,cis,xcs,ycs,rfulls,rrs,t,chis,chitrs] = ensemble_del_psi(savename,Nits,discardfactor,N,dt,steps,L,plotevery,tau,svar,kappa,lzero,tavg,Dthetas,Sgrad)
% function [vxs,vys,cis,xcs,ycs,rfulls,rrs,t,chis,chitrs] = ensemble_del_psi(savename,Nits,discardfactor,N,dt,steps,L,plotevery,tau,svar,kappa,lzero,tavg,Dthetas,Sgrad)


if(isdeployed)

   Nits = str2num(Nits)
   discardfactor = str2num(discardfactor)
   N = str2num(N)
   dt = str2num(dt)
   steps = str2num(steps)
   %vr = str2num(vr)
   L = str2num(L)
   plotevery = str2num(plotevery)
   tau = str2num(tau)
   svar = str2num(svar)
   kappa = str2num(kappa)
   lzero = str2num(lzero)
   tavg = str2num(tavg)
   Dthetas = str2num(Dthetas)
   Sgrad = str2num(Sgrad)    

end

rng('shuffle');
scurr = rng


dodraw = 0;


saveeveryits = round(100*(steps/2e4))


vxs = zeros(length(Dthetas),Nits);
vys = zeros(length(Dthetas),Nits);
cis = zeros(length(Dthetas),Nits);
cits = zeros(length(Dthetas),Nits);
xcs = cell(length(Dthetas),Nits);
ycs = cell(length(Dthetas),Nits);
rlams = cell(length(Dthetas),Nits);
rfulls = cell(length(Dthetas),Nits);
rfullsxx = cell(length(Dthetas),Nits);
rfullsyy = cell(length(Dthetas),Nits);
rrs = cell(length(Dthetas),Nits);
chis = NaN*ones(1,Nits);
chitrs = NaN*ones(1,Nits);
%Amats = zeros(3,Nits); % Axx, Axy, Ayy
Aimall = zeros(2,2,Nits);
gxvar = zeros(1,Nits);
gyvar = zeros(1,Nits); % this will be nuts if length(Dthetas) not 1
gxTvar = zeros(1,Nits);
gyTvar = zeros(1,Nits); % this will be nuts if length(Dthetas) not 1

if(length(Dthetas)>1)
  warning('Probably too many dthetas, a problem')
end


%discardsteps = ceil(tavg*discardfactor/(dt*plotevery));
discardsteps = ceil(max(tavg,tau)*discardfactor/(dt*plotevery)); % drop off discartfactor x either the correlation time or the averaging time

gxss = zeros(Nits,steps/plotevery-discardsteps);
gyss = zeros(Nits,steps/plotevery-discardsteps);
gxTss = zeros(Nits,steps/plotevery-discardsteps);
gyTss = zeros(Nits,steps/plotevery-discardsteps);


for i = 1:length(Dthetas)
    
    for s = 1:Nits
        try
            [gxs,gys,xss,yss,vxss,vyss,thetass,phis,systematic,gxTs,gyTs] = rotor_del_psi(N,dt,steps,L,plotevery,dodraw,tau,svar,kappa,lzero,tavg,Dthetas(i),Sgrad);
        %[gxs,gys,xss,yss,vxss,vyss,thetass,phiss,systematic] = rotor_align_estimate(N,dt,steps,vr,L,plotevery,dodraw,tau,svar,chi,ell,tavg,Dthetas(i),Sgrad);
        catch err        
            getReport(err)
            [gxs,gys,xss,yss,vxss,vyss,thetass,phis,systematic] = rotor_del_psi(N,dt,steps,L,plotevery,dodraw,tau,svar,kappa,lzero,tavg,Dthetas(i),Sgrad);
        end
        xss = xss(:,discardsteps+1:end);
        yss = yss(:,discardsteps+1:end);
        gxs = gxs(:,discardsteps+1:end);
        gys = gys(:,discardsteps+1:end);
	gxTs = gxTs(:,discardsteps+1:end);
        gyTs = gyTs(:,discardsteps+1:end);

gxss(s,:) = gxs;
gyss(s,:) = gys;
gxTss(s,:) = gxTs;
gyTss(s,:) = gyTs;

Aim = Ainv_means_pos(xss,yss);
Aimall(:,:,s) = Aim;

gxvar(s) = mean((gxs(:)-Sgrad).^2);
gyvar(s) = mean((gys(:)).^2);
gxTvar(s) = mean((gxTs(:)-Sgrad).^2);
gyTvar(s) = mean((gyTs(:)).^2);
        cutpoint = 1;
        maxlag = length(xss)-1;
        donorm = 0; % we want to collect the *un-normalized* lambda correlations, so when we average them we are averaging sanely
        %[rlam,rr,t] = g_and_r_corrs_updated(gxs-Sgrad,gys,xss,yss,cutpoint,maxlag,dt,plotevery,donorm);
%[t,rr,rrs_inst,rf,rfs,chi_inst,chi_tr_inst] = full_corr(xss,yss,cutpoint,maxlag,dt,plotevery);
%[t,rr,rrs
[t,rr,rrs_inst,rf,rfs,chi_inst,chi_tr_inst,rfxx,rfyy,rfxxs,rfyys] = full_corr_xy(xss,yss,cutpoint,maxlag,dt,plotevery);
        rfulls{i,s} = rf;
	rfullsxx{i,s} = rfxx;
	rfullsyy{i,s} = rfyy;
        rrs{i,s} = rr;
        vxss = vxss(:,discardsteps+1:end);
        vyss = vyss(:,discardsteps+1:end);
        

chis(s) = mean(chi_inst);
chitrs(s) = mean(chi_tr_inst);
        %chis(s) = 0.5*mean(sum(dxss.^2+dyss.^2));
	%chitrs(s) = mean(chi_trace(dxss,dyss));


        
        [ci,cit] = cilength(xss,yss);
        cis(i,s) = ci;
        cits(i,s) = cit;
        vxs(i,s) = mean(vxss(:));
        vys(i,s) = mean(vyss(:));
        xcs{i,s} = mean(xss)-mean(xss(:,1));
        ycs{i,s} = mean(yss)-mean(yss(:,1));
        fprintf('%d ',s);
   if(rem(s,saveeveryits)==0)
     save(savename)
   end
    end
    fprintf('\n');
    fprintf('%d its w/Dth = %3.2g, ci = %3.2g +/- %3.2g, vx = %3.2g +/- %3.2g \n',Nits,Dthetas(i),mean(cis(i,:)),std(cis(i,:))/sqrt(Nits),mean(vxs(i,:)),std(vxs(i,:))/sqrt(Nits));
    save(savename)
end
