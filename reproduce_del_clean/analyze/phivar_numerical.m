function [dphis,cosm] = phivar_numerical(sigmas)

dphis = zeros(size(sigmas));
cosm = zeros(size(sigmas));

for s = 1:length(sigmas(:))
try
sig = sigmas(s);
    %fi = @(x,y) exp(-(x.^2)/(2*sig^2)).*exp(-(y.^2)/(2*sig^2));
    fi = @(x,y) exp(-(x.^2)/(2*sig^2)).*exp(-(y.^2)/(2*sig^2)).*(atan(y./(1+x)).^2);    
    dphivar1 = integral2(fi,-inf,inf,-inf,inf);
    dphivar2 = integral2(fi,-4*sig,4*sig,-4*sig,4*sig);
    dphivar = max(dphivar1,dphivar2);  % ok since we know that the integrand is positive everywhere - in principle should always be the 
                                       % top one, but the adaptive
                                       % quadrature can screw up
   
    dphivar = dphivar/(2*pi*sig^2);
    
    dphis(s) = sqrt(dphivar);
    fc = @(x,y) exp(-(x.^2)/(2*sig^2)).*exp(-(y.^2)/(2*sig^2)).*(1+x)./sqrt((1+x).^2+y.^2);                                     
    cm1 = integral2(fc,-inf,inf,-inf,inf);
    cm1 = cm1/(2*pi*sig^2);
    cosm(s) = cm1;
    %dphivar2 = integral2(fi,-4*sig,4*sig,-4*sig,4*sig,'RelTol',1e-10,'AbsTol',1e-10);
   catch err
   getReport(err)
   cosm(s) = NaN;
   dphis(s) = NaN;
   end
end