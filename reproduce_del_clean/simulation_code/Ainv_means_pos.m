function [Aim,Aims] = Ainv_means_pos(xss,yss)

  Aims = zeros(2,2,size(xss,2));

for k = 1:size(xss,2)
	  xs = xss(:,k); ys = yss(:,k);
dxs = xs-mean(xs); dys = ys-mean(ys);

    
    A = [sum(dxs.^2) sum(dxs.*dys) ; sum(dxs.*dys) sum(dys.*dys)];
    Ainv = inv(A);
    Aims(:,:,k) = Ainv;

end

Aim = mean(Aims,3);

end
