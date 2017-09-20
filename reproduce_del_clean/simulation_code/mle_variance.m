function gv = mle_variance(xs,ys)
%% assumes isotropy version:
% A = mob_matrix_norm(xs,ys,Rz);
% chi = 0.5*trace(A);
% [qx,qy] = xdirs(xs,ys,Rz);
% N = length(xs);
% gv = (1/chi^2)*sum(qx.^2+qy.^2)/(N^2);
xs = xs-mean(xs); ys = ys-mean(ys);
Amat = [sum(xs.^2) sum(xs.*ys) ; sum(ys.*xs) sum(ys.^2)];
Ainv = inv(Amat);
% qx = xs-mean(xs); qy = ys-mean(ys);
% B = [sum(qx.^2) sum(qx.*qy) ; sum(qx.*qy) sum(qy.^2)];
% C = Ainv*B;
%gv = sum(sum(Ainv.*C));
gv = trace(Ainv);
%ghat = linsolve(A,z);
%gv = ghat(1)^2+ghat(2)^2;