function chitr = chi_trace(xss,yss)

chitr = zeros(1,length(xss));

for k = 1:size(xss,2)
    gv = mle_variance(xss(:,k),yss(:,k));
    chitr(k) = 2/gv;

end