function c = cos_average(sigmas)
c = zeros(size(sigmas));
for s = 1:length(sigmas)
    sig = sigmas(s);
    
    if((sig<10)&&(sig>0.05))
        fc = @(x,y) exp(-(x.^2)/(2*sig^2)).*exp(-(y.^2)/(2*sig^2)).*(1+x)./sqrt((1+x).^2+y.^2);                                     
        cm1 = integral2(fc,-inf,inf,-inf,inf,'RelTol',1e-10,'AbsTol',1e-10);
        cm1 = cm1/(2*pi*sig^2);
        c(s) = cm1;
    elseif(sig>=1e1)
        c(s) = sqrt(pi/8)/sig; % asymptotic value
    elseif(sig<=0.05)
        c(s) = 1-(sig^2)/2;
    else
        fprintf('May have a problem with sig = %3.6g',sig);
        break
    end
    %dphivar2 = integral2(fi,-4*sig,4*sig,-4*sig,4*sig,'RelTol',1e-10,'AbsTol',1e-10);
    
end
end
