function [xs,ys] = hex_remove_start(N)

b = 3; a = 3;
c = (1-N);
q1 = (-b + sqrt(b^2-4*a*c))/(2*a);
q2 = (-b - sqrt(b^2-4*a*c))/(2*a);
qm = max(q1,q2);
q = ceil(qm);
[xs,ys] = hex_packed(1,q);
nr = length(xs)-N;
bad = randperm(length(xs),nr);
xs(bad) = [];
ys(bad) = [];
