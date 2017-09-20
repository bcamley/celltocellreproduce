function LL = loglike_mult(gx,gy,M,xs,ys,nr,KD,sigmad)

mu = 1 + gx*xs + gy*ys;
ci = max(mu,0); % we ensure positivity here so nothing blows up too much?
deltac_squared = (1/nr)*((ci/KD+1).^2).*(ci*KD);
h = deltac_squared + sigmad^2;

LL = -0.5*sum(log(h))-0.5*sum( ((mu-M).^2)./h);

end