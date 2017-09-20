function [gv,gxm,gym,gall] = full_deltac_deltasig(Nits,Ns,gx,gy,nr,KDs,sigmads)
% function [gv,gxm,gym] = full_deltac_deltasig(Nits,Ns,gx,gy,nr,KDs,sigmads)

gv = zeros(length(Ns),length(sigmads),length(KDs));
gxm = zeros(length(Ns),length(sigmads),length(KDs));
gym = zeros(length(Ns),length(sigmads),length(KDs));
gall = cell(length(Ns),length(sigmads),length(KDs));

for i = 1:length(Ns)
    
    for j = 1:length(sigmads)
    	for k = 1:length(KDs)
	        [xs,ys] = hex_remove_start(Ns(i));
        	[gxest,gyest] = mult_test_conc(Nits,xs,ys,gx,gy,nr,KDs(k),sigmads(j));
		        gv(i,j,k) =  mean((gxest-gx).^2+(gyest-gy).^2);
        		gxm(i,j,k) = mean(gxest);
        		gym(i,j,k) = mean(gyest);
			gall{i,j,k} = [gxest ; gyest];
        
    end
    end
    fprintf('Ns %d /%d, sds %d / %d, KDs = %d/%d, gv = %3.6g \n',i,length(Ns),j,length(sigmads),k,length(KDs),gv(i,j,k));
end

end

