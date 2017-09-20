function [gxest,gyest,Ms,ggs] = mult_test_conc(Nits,xs,ys,gx,gy,nr,KD,sigmad)
% function [gxest,gyest] = mult_test_conc(Nits,xs,ys,gx,gy,nr,KD,sigmad)
gxest = zeros(1,Nits);
gyest = zeros(1,Nits);
ggs = cell(1,Nits);
Ms = cell(1,Nits);
%opts = optimset('Display','off','TolX',1e-3);
mu = 1 + gx*xs + gy*ys;
ci = mu;
deltac_squared = (1/nr)*((ci/KD+1).^2).*(ci*KD);
deltac_err = sqrt(deltac_squared);

for i = 1:Nits
    M = (1+gx*xs+gy*ys) + deltac_err.*randn(size(xs)) + sigmad*randn(size(xs));
    Ms{i} = M;
    [ghats,ghat_guess] = ghat_conc_and_sigmad(xs,ys,M,nr,KD,sigmad);
    gxest(i) = ghats(1);
    gyest(i) = ghats(2);
    ggs{i} = ghat_guess;
end

end
