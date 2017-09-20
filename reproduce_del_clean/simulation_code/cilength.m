function [ci,cit] = cilength(xss,yss)

if(size(xss,1)>1)
    xs = mean(xss);
    ys = mean(yss);
else
    xs = xss;
    ys = yss;
end

dx = diff(xs);
dy = diff(ys);
ds = sqrt(dx.^2+dy.^2);

ci = trapz(dx)/trapz(ds);
cit = mean(dx./ds);