function od = cellvar_analyze(filename_start,Nboots,doplot,doCIerr,nw,nh)
% od = cellvar_analyze(filename_start,Nboots,doplot,doCIerr,nw,nh)
%
% filename_start is the first part of the files to be analyzed, i.e. filename_start1.mat, filename_start2.mat, ...
% Nboots is the number of bootstrap runs used in estimating the errors in the position-position correlation time (default 50)
% if doplot is true, after finishing, the simulation will plot position-position correlations and their exponential fits
% doCIerr should be false - it is an attempt to propagate some errors through the CI, but turns out to be both slow and negligible
% nw and nh are the number of points simulated for each parameter (i.e. the numbers given to twosweep.pl)
%

od = struct;
nmax = nw*nh;

citsall = NaN*ones(1,nmax); % time-averaged cluster chemotactic index
citseall = NaN*ones(1,nmax); % error estimate on cit
Dtheta = citsall;	     % parameter Dtheta for simulation
tavgs = citsall;	     % parameter Tavg (averaging time) for simulation
Nitcount = citsall;	     % number of iterations done for each simulation
trrs = citsall;		     % position-position correlation times
trr_errs = citsall;	     % bootstrap error in position-position correlation times
trr_fulls = citsall;	     % 
trr_err_fulls = citsall;     % 	 
ts = cell(1,nmax);	     % time window for delta r autocorrelation
rras = cell(1,nmax);	     % delta r autocorrelation
Ns = citsall;		     % number of cells for simulations
kappas = citsall;	     % stiffness/adhesion kappa for each simulation
chisall = citsall;	     % chi for the cell configuration
chisall_err = citsall;	     % estimate of statistical error in chi
chitrsall = citsall; 

gxvs = citsall;
gyvs = citsall;
gxTvs = citsall;
gyTvs = citsall;
Aim = cell(1,nmax);

%nn = round(sqrt(nmax));

od.xss = cell(1,nmax);
od.yss = cell(1,nmax);
od.thetass = cell(1,nmax);
od.systematic = cell(1,nmax);

if(doplot)
	clf;

end

for s = 1:nmax
   filename = sprintf('%s%d.mat',filename_start,s);
   try
   load(filename,'cis','Dthetas','tavg','Nits','rlams','rrs','t','N','cits','chis','svar','Sgrad','kappa','rfulls','chitrs','gxvar','gyvar','gxTvar','gyTvar','Aimall','xss','yss','thetass','systematic','Sgrad');
   od.Sgrad = Sgrad; % should be constant...
   od.xss{s} = xss(:,end);
   od.yss{s} = yss(:,end);
   od.thetass{s} = thetass(:,end);
   od.systematic{s} = systematic;
   chisall(s) = mean(chis);
       chitrsall(s) = mean(chitrs);
       chisall_err(s) = std(chis)/sqrt(Nits);
       kappas(s) = kappa;
       Ns(s) = N;

       Aim{s} = mean(Aimall,3);

       gxvs(s) = mean(gxvar);
       gyvs(s) = mean(gyvar);
       gxTvs(s) = mean(gxTvar);
       gyTvs(s) = mean(gyTvar);

       
       citsall(s) = mean(cits); % will probably throw an error if we ran more than one Dtheta per site
       Dtheta(s) = Dthetas;
       tavgs(s) = tavg;  
       Nitcount(s) = Nits;
       citseall(s) = std(cits)/sqrt(Nitcount(s)); % this is only a very rough estimate, of course... 
	 told = t;
       [tau,tau_err,taus,t,rra_true] = autocorr_fit_bootstrap(rrs,told,Nboots);  
       trrs(s) = tau;
       trr_errs(s) = tau_err;
       rras{s} = rra_true;
       ts{s} = t;
       [tau_full,tau_err_full,taus_full,t,rra_true_full] = autocorr_fit_bootstrap(rfulls,told,Nboots);  
       trr_fulls(s) = tau_full;
       trr_err_fulls(s) = tau_err_full;
rrfs{s} = rra_true_full;

%        try
%                  tau_fullx = autocorr_fit_bootstrap(rfullsxx,told,Nboots);  
%                  tau_fully = autocorr_fit_bootstrap(rfullsyy,told,Nboots);  
% od.trrx(s) = tau_fullx;
% od.trry(s) = tau_fully;
%        catch err
%        getReport(err)
%        disp('Problem with x/y correlations')
%        end

       fprintf('(%d/%d): Dtheta %3.3g, tavg %3.3g, kappa %3.3g, N %d, ci = %3.3g +/- %3.3g \n',s,nmax,Dthetas,tavg,kappa,N,citsall(s),citseall(s));
              if(doplot)
            subplot(nw,nh,s);
            plot(t,rra_true_full,'k');
            hold on
            plot(t,exp(-t/trr_fulls(s)),'b--');
        end
   catch err
       trrs(s) = Inf;	 
       trr_fulls(s) = Inf;	 
       fprintf('Error on analyzing file %s \n',filename);
       getReport(err)
   end
    
end


try
Dtheta = reshape(Dtheta,nw,nh);
tavgs = reshape(tavgs,nw,nh);
Nitcount = reshape(Nitcount,nw,nh);
trrs = reshape(trrs,nw,nh);
trr_errs = reshape(trr_errs,nw,nh);
chisall = reshape(chisall,nw,nh);
chitrsall = reshape(chitrsall,nw,nh);
chisall_err = reshape(chisall_err,nw,nh);
citsall = reshape(citsall,nw,nh);
citseall = reshape(citseall,nw,nh);
trr_fulls = reshape(trr_fulls,nw,nh);
trr_err_fulls = reshape(trr_err_fulls,nw,nh);
gxTvs = reshape(gxTvs,nw,nh);
gyTvs = reshape(gyTvs,nw,nh);
gxvs = reshape(gxvs,nw,nh);
gyvs = reshape(gyvs,nw,nh);

kappas = reshape(kappas,nw,nh);
Ns = reshape(Ns,nw,nh);

catch err
    getReport(err)
end

try
SNR = (Sgrad^2)./(svar^2./chisall);
SNRT = SNR.*(1+tavgs./trrs);

SNRfull = (Sgrad^2)./(svar^2./chitrsall);
SNRTfull = SNRfull.*(1+tavgs./trr_fulls);

disp('Evaluating predicted CI (slowish)...')

[dphis,cosm] = phivar_numerical(1./sqrt(SNRT));
[dphisfull,cosmfull] = phivar_numerical(1./sqrt(SNRTfull));

if(doCIerr)
	%
	mco = zeros(size(SNRT));
	msco =  zeros(size(SNRT));


	fprintf('Bootstrapping CIs (VERY slow), b = ')
	for b = 1:Nboots
		SNRb = (Sgrad^2)./(svar^2./(chisall+randn(size(chisall)).*chisall_err));
		SNRTb = SNRb.*(1+tavgs./(trrs+randn(size(trrs)).*trr_errs));
		[dphisb,cosmb] = phivar_numerical(1./sqrt(SNRTb));
		mco = mco + cosmb;
		msco = msco + cosmb.^2;
		fprintf('%d ',b);
	end
	fprintf('\n');
	mco = mco/Nboots;
	msco = msco/Nboots;
	od.cosm_std = sqrt(msco-mco.^2);

end

catch err
    getReport(err)
    cosm = NaN;
end

%[Dtheta,tavgs,Nitcount,trrs,chisall,chisall_err,citsall,citseall,SNR,SNRT,cosm,ts,rras,Ns,kappas] = cellvar_cleaned_bootstrap(filename_start,nmax,doplot,doCIerr)
od.Dtheta = Dtheta;
od.tavgs = tavgs;
od.Nitcount = Nitcount;
od.trrs = trrs;
od.chisall = chisall;
od.chisall_err = chisall_err;
od.citsall = citsall;
od.citseall = citseall;
od.SNR = SNR;
od.SNRT = SNRT;
od.cosm = cosm;
od.SNRfull = SNRfull;
od.SNRTfull=SNRTfull;
od.cosmfull = cosmfull;
od.ts = ts;
od.rras = rras;
od.Ns = Ns;
od.kappas = kappas;
od.trr_errs = trr_errs;
od.trr_fulls = trr_fulls;
od.trr_err_fulls = trr_err_fulls;
od.chitrsall = chitrsall;
od.rrfs = rrfs;
od.gxvs = gxvs;
od.gyvs = gyvs;
od.gxTvs = gxTvs;
od.gyTvs = gyTvs;
od.Aim = Aim;
