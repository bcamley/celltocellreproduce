function [t,rr,rrs,rf,rfs,chis,chi_trs,rfxx,rfyy,rfxxs,rfyys] = full_corr_xy(xss,yss,cutpoint,maxlag,dt,plotevery)

xss = xss(:,cutpoint:end);
yss = yss(:,cutpoint:end);

dxss = xss-ones(size(xss,1),1)*mean(xss);
dyss = yss-ones(size(yss,1),1)*mean(yss);

Axx = sum(dxss.*dxss);
Axy = sum(dyss.*dxss);
Ayy = sum(dyss.*dyss);

rrs = zeros(size(xss,1),2*maxlag+1);
rfs = zeros(size(xss,1),2*maxlag+1);
rfxxs = zeros(size(xss,1),2*maxlag+1);
rfyys = zeros(size(xss,1),2*maxlag+1);

lamaxs = zeros(size(dxss)); % averaged lambda-ish vector
lamays = zeros(size(dxss)); % averaged lambda-ish vector

chis = 0.5*(Axx+Ayy);
chi_trs = zeros(size(chis));

for s = 1:size(xss,2)
	  A = [Axx(s) Axy(s) ; Axy(s) Ayy(s)];
chi_trs(s) = 2/trace(inv(A));
	  for n = 1:size(xss,1)
		    lama = linsolve(A,[dxss(n,s) ; dyss(n,s)]);
lamaxs(n,s) = lama(1);
lamays(n,s) = lama(2);
          end
end
		    

for s = 1:size(xss,1)
    rs = xcorr([dxss(s,:) ; dyss(s,:)].',maxlag,'unbiased');
    lamcorr = xcorr([lamaxs(s,:) ; lamays(s,:)].',maxlag,'unbiased');
    rsum = (rs(:,1)+rs(:,4));
    rrs(s,:) = rsum;
    rfsum = (lamcorr(:,1)+lamcorr(:,4));
    rfs(s,:) = rfsum;
    rfxxs(s,:) = lamcorr(:,1);
    rfyys(s,:) = lamcorr(:,4);
end


rr = mean(rrs);
rf = mean(rfs);
rfxx = mean(rfxxs);
rfyy = mean(rfyys);

t = dt*plotevery*(-maxlag:maxlag);
