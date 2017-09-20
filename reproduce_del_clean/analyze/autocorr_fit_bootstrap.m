function [tau,tau_err,taus,t,rra_true] = autocorr_fit_bootstrap(rrs,t,Nboots)



    rra = mean(cell2mat(rrs.'));
    maxlag = (length(rra)-1)/2;
    t = t(maxlag+1:end);
        rra = rra(maxlag+1:end);
       rra = rra/rra(1);
    rra_true = rra;

taus = zeros(1,Nboots);
for b = 1:Nboots
boot = randi(length(rrs),[1 length(rrs)]); % choose random sample of iterations, with replacement 
    rrb = rrs(boot);
    rra = mean(cell2mat(rrb.'));
       rra = rra(maxlag+1:end);
       rra = rra/rra(1);
       kr = find(rra<0.8,1);
       kr = min(kr*4,length(t)); % a rough fitting window
       tf = t(1:kr);
       taus(b) = nlinfit(tf,rra(1:kr),@(beta,x) exp(-x/beta),t(kr)*2);
end

tau = mean(taus);
tau_err = std(taus);

end
