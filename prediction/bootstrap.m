rng default
%bootstrap estimation
% mean_xbs=[];std_xbs=[];
parfor sim=1:300
% sim
    xi_xb=xi_xb_est{sim};
    xi_xs=xi_xs_est{sim};
    for i=1:2
        for j=1:2
            [bootstat,bootsam] = bootstrp(1000,@cov,xi_xb(:,i),xi_xs(:,j));
            mean_tmp=mean(bootstat);
            std_tmp=std(bootstat);
            mean_xbs(i,j,sim)=mean_tmp(2);
            std_xbs(i,j,sim)=std_tmp(2);
        end
    end
end
mean(mean_xbs,3)
mean(std_xbs,3)

parfor sim=1:300

    xi_yb=xi_yb_est{sim};
    xi_ys=xi_ys_est{sim};
    for i=1:2
        for j=1:2
            [bootstat,bootsam] = bootstrp(1000,@cov,xi_yb(:,i),xi_ys(:,j));
            mean_tmp=mean(bootstat);
            std_tmp=std(bootstat);
            mean_ybs(i,j,sim)=mean_tmp(2);
            std_ybs(i,j,sim)=std_tmp(2);
        end
    end
end
mean(mean_ybs,3)
mean(std_ybs,3)

%true estimation
for i=1:rep
%     i
    xi_xb_tmp=xi_xb_est{i};
    xi_xs_tmp=xi_xs_est{i};
    cov_xix(:,:,i)=cov([xi_xb_tmp(:,1:2) xi_xs_tmp(:,1:2)]);

    xi_yb_tmp=xi_yb_est{i};
    xi_ys_tmp=xi_ys_est{i};
    cov_xiy(:,:,i)=cov([xi_yb_tmp(:,1:2) xi_ys_tmp(:,1:2)]);
   
end

cov_xix_mean=mean(cov_xix,3);
cov_xiy_mean=mean(cov_xiy,3);
cov_xix_std=[std(cov_xix(1,3,:)) std(cov_xix(1,4,:)); std(cov_xix(2,3,:)) std(cov_xix(2,4,:))];
cov_xiy_std=[std(cov_xiy(1,3,:)) std(cov_xiy(1,4,:)); std(cov_xiy(2,3,:)) std(cov_xiy(2,4,:))];
[cov_xix_mean(1:2,3:4) cov_xiy_mean(1:2,3:4)]
[cov_xix_std cov_xiy_std]



