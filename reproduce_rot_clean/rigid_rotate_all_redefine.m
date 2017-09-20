function [V,om,SS,TT,c_interp,sigma_interp] = rigid_rotate_all_redefine(SNR,tau,ninterp,np)
% function [V,om,SS,TT] = rigid_rotate_all_redefine(SNR,tau,ninterp,np)
% here SNR is SNR at zero rotation = (g^2)*(chi)/(sigma_delta^2)
%                  where g is the gradient, chi = 0.5*sum_i (x_i^2+y_i^2)
%                  (assuming isotropy), sigma_delta the std of the
%                  cell-cell variance
%
% tau = (averaging time)/(rotation time) = (averaging time)/(R/vmax)
% 

[SS,TT] = meshgrid(SNR,tau);
V = NaN*ones(size(SS));
om = NaN*ones(size(SS));

sig_max = 2*sqrt(1/min(SNR))
sigma_interp = linspace(0.1/np,sig_max,ninterp);
c_interp = cos_average(sigma_interp);

for i = 1:length(SNR)
    for j = 1:length(tau)
        om_max = sqrt(2);
        oms = linspace(0,om_max,np);
        sigmas = sqrt((1/SNR(i))*1./(1+(tau(j)*oms).^2));
        c = cos_average_interp(sigmas,sigma_interp,c_interp);
        vmags = sqrt(1-0.5*(oms.^2));
        vs = c.*vmags;
        [vmax,k] = max(vs);
        V(j,i) = vmax;
        om(j,i) = oms(k);
        
    end
    

end

end

