clear all
load data300_mn200_bs
load DVb_mn200
load DVs_mn200
load cov_xi_xy_mn200

rep=300;
n=200;m=n;
grid_length=100;

%% estimating equation for gamma of cov(xi_zb,xi_zs)
for i=1:rep
n_covz=2;
    sigma0=ones(1,n_covz^2)*0.1;
    data_b=DATA_B{i};
    data_s=DATA_S{i};
    
    [pair_ij,f2]=pair_eez_fun(data_b,data_s,Vz_b(:,:,i),Vz_s(:,:,i),n_covz,n,m,grid_length);
    options = optimoptions('fsolve','Display','off');
    [paraz(i,:),fval,exitflag(i,:),output] = fsolve(@eez_fun,sigma0,options,n,m,n_covz,f2,pair_ij,cov_xix(:,:,i),cov_xiy(:,:,i),...
        Vx_b(:,:,i),Vy_b(:,:,i),Vz_b(:,:,i),Vx_s(:,:,i),Vy_s(:,:,i),Vz_s(:,:,i));
    
end
save cov_xiz_mn200 paraz

%% estimate cov_xix, cov_xiy,
load xix_logis_sim300_mn400
load xiy_logis_sim300_mn400

for i=1:rep
%     i
    xi_xb_tmp=xi_xb_est{i};
    xi_xs_tmp=xi_xs_est{i};
    cov_xix(:,:,i)=cov([xi_xb_tmp(:,1:2) xi_xs_tmp(:,1:2)]);

    xi_yb_tmp=xi_yb_est{i};
    xi_ys_tmp=xi_ys_est{i};
    cov_xiy(:,:,i)=cov([xi_yb_tmp(:,1:2) xi_ys_tmp(:,1:2)]);
   
end
save cov_xi_xy_mn400 cov_xix cov_xiy

load cov_xi_xy_mn400
cov_xix_mean=mean(cov_xix,3);
cov_xiy_mean=mean(cov_xiy,3);
cov_xix_std=[std(cov_xix(1,3,:)) std(cov_xix(1,4,:)); std(cov_xix(2,3,:)) std(cov_xix(2,4,:))];
cov_xiy_std=[std(cov_xiy(1,3,:)) std(cov_xiy(1,4,:)); std(cov_xiy(2,3,:)) std(cov_xiy(2,4,:))];
[cov_xix_mean(1:2,3:4) cov_xiy_mean(1:2,3:4)]
[cov_xix_std cov_xiy_std]

%combine two files
% load xix1_logis_sim150_mn400 
% load xix2_logis_sim150_mn400 
% for i=1:150
% xi_xb_all{i}=xi_xb_est{i};
% xi_xs_all{i}=xi_xs_est{i};
% end
% xi_xb_est=xi_xb_all;
% xi_xs_est=xi_xs_all;
% save xix_logis_sim300_mn400 xi_xb_est xi_xs_est
% 
% load xiy1_logis_sim150_mn400 
% load xiy2_logis_sim150_mn400 
% for i=1:150
% xi_yb_all{i}=xi_yb_est{i};
% xi_ys_all{i}=xi_ys_est{i};
% end
% xi_yb_est=xi_yb_all;
% xi_ys_est=xi_ys_all;
% save xiy_logis_sim300_mn400 xi_yb_est xi_ys_est
% 
% 
