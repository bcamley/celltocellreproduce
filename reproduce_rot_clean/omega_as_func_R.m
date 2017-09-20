function OMEGA = omega_as_func_R(Rs,T,sigmad,rhoc,g,vmax,sigma_interp,c_interp,np)
  chi = (pi/4)*rhoc*Rs.^4;
SNR = (g^2)*chi/(sigmad^2);
omega = zeros(size(Rs));
tau = T*vmax./Rs;
for i = 1:length(SNR)
        om_max = sqrt(2);
        oms = linspace(0,om_max,np);
        sigmas = sqrt((1/SNR(i))*1./(1+(tau(i)*oms).^2));
        c = cos_average_interp(sigmas,sigma_interp,c_interp);
        vmags = sqrt(1-0.5*(oms.^2));
        vs = c.*vmags;
        [vmax,k] = max(vs);
omega(i) = oms(k);
end

OMEGA = omega*vmax./Rs; % the unitful capital Omega

end
