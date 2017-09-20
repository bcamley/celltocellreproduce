function c = cos_average_interp(sigmas,sigma_interp,c_interp)
    c = interp1(sigma_interp,c_interp,sigmas);
    upasy = sigmas>=10;
    lowasy = sigmas<=0.05;
    c(upasy) = sqrt(pi/8)./sigmas(upasy);
    c(lowasy) = 1-(sigmas(lowasy).^2/2);
    bad = isnan(c);
    if(sum(bad)>0)
      fprintf('Problem with interpolation, rerunning %d points, could be very slow \n',length(sigmas(bad)))
       %disp(sigmas(bad))
       c_new = cos_average(sigmas(bad));
       c(bad) = c_new;
    end

    

end


