function [rg,rr,t,rrs,rmsq] = g_and_r_corrs_updated(gxs,gys,xss,yss,cutpoint,maxlag,dt,plotevery,donorm)

%gxs = gxs(1:plotevery:end);
%gys = gys(1:plotevery:end);
gxs = gxs(cutpoint:end);
gys = gys(cutpoint:end);

xss = xss(:,cutpoint:end);
yss = yss(:,cutpoint:end);

dxss = xss-ones(size(xss,1),1)*mean(xss);
dyss = yss-ones(size(yss,1),1)*mean(yss);

rmsq = mean(mean(dxss.^2+dyss.^2));

rgun = xcorr([gxs ; gys].',maxlag,'unbiased');
if(donorm)
    rg = (rgun(:,1)+rgun(:,4))/(rgun(maxlag+1,1)+rgun(maxlag+1,4));  % normalized
else
    rg = rgun(:,1)+rgun(:,4);
end

rrs = zeros(size(xss,1),2*maxlag+1);

for s = 1:size(xss,1)
    rs = xcorr([dxss(s,:) ; dyss(s,:)].',maxlag,'unbiased');
    rsum = (rs(:,1)+rs(:,4));
    rrs(s,:) = rsum;
end

%size(rrs)
rrun = mean(rrs);
if(donorm)
    rr = rrun/rrun(maxlag+1);  % normalized
else
    rr = rrun;
end
t = dt*plotevery*(-maxlag:maxlag);