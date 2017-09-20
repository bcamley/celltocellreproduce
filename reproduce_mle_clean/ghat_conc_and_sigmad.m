function [ghats,ghat_guess] = ghat_conc_and_sigmad(xs,ys,M,nr,KD,sigmad)
  % note that this assumes xs,ys is centered!!
%ghat_guess = ghat_multiplicative(xs,ys,M,a,b);
Amat = [ sum(xs.^2) sum(xs.*ys) ; sum(xs.*ys) sum(ys.^2)];
ghat_guess = linsolve(Amat,[sum(M.*xs) ; sum(M.*ys)]);

% a less informative initial guess will also work
%LLg = -Inf; lc = 0;
%while(isinf(LLg))%
%	ghat_guess = [0.2*randn 0.2*randn]; % localized guess
%	LLg = -loglike_mult(ghat_guess(1),ghat_guess(2),M,xs,ys,nr,KD,sigmad);
%	lc = lc + 1;
%	if(lc>10)
%		error('inflike') % can't find a non-infinite likelihood to star%t out... can happen if initial guesses are very bad
%	end
%end

opts = optimset('Display','off','TolX',1e-3,'MaxFunEvals',5e4,'MaxIter',5e4);
ghats = fminsearch(@(g) -loglike_mult(g(1),g(2),M,xs,ys,nr,KD,sigmad),ghat_guess,opts);


end


