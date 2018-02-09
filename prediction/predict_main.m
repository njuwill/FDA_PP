%% logistic likelihood method
clear all
load('C:\Users\mwang\Dropbox\finance project\pegasus\data\sim300\data300_nm100_rho5.mat')
load cov_sepH_sim300_mn100
rep=100;n=100;m=n;grid_length=100;
for i=1:rep
   i
   data_b=DATA_B{i};
   data_s=DATA_S{i};
   [xi_xb_est{i},flagx]=predict_logisticx_fun(data_b,n,grid_length,cov_xb(:,:,i));
   [xi_xs_est{i},flagx]=predict_logisticx_fun(data_s,n,grid_length,cov_xs(:,:,i));
end
save xix_logis_sim300_mn200 xi_xb_est xi_xs_est

% predict xi_y without considering day lags
for i=1:rep
   data_b=DATA_B{i};
   data_s=DATA_S{i};
   [xi_yb_est{i}]=predict_logisticy_fun(data_b,m,grid_length,cov_yb(:,:,i));
   [xi_ys_est{i}]=predict_logisticy_fun(data_s,m,grid_length,cov_ys(:,:,i));
end
save xiy_logis_sim300_mn200 xi_yb_est xi_ys_est

% predict xi_y with considering day lags
for i=1:rep
   i
   data=DATA{i};
   [paray_logistic_lag{i},flagy{i}]=predict_logisticy_lag_fun(data,m,grid_length,cov_y(:,:,i));
end
save xiy_logis_lag_sim100_mn100 paray_logistic

%% maximize likelihood method mn100 case
load('C:\Users\mwang\Dropbox\finance project\pegasus\data\sim300\data300_nm100_rho5.mat')
% load data100_nm200_rho5.mat
grid_length=100;
T=1;
rep=300;
m=100;n=m;
tic
load('C:\Users\mwang\Dropbox\finance project\matlab mfpca\simplified h selection\h with z2 direction\h100')
load('C:\Users\mwang\Dropbox\finance project\matlab mfpca\simplified decomposition\result sim300 z2 directions\cov_sepH_sim300_mn100')
for i=1:rep
   i
   data=DATA{i};
   hrho=hrho_100(i);
   [parax,flagx,paray,flagy]=predict_xy_fun(data,m,n,grid_length,hrho,cov_x(:,:,i),cov_y(:,:,i));
end
toc
[mean(parax);var(parax)]
[mean(paray);var(paray)]

save xix_sim300_mn100 parax flagx


%% new estimating equation method, without estimating first order intensity
for i=1:rep
   i
   data=DATA{i};
   [parax_ee{i},flagx]=predict_x_fun(data,n,grid_length,cov_x(:,:,i));
end
save xix_ee_sim300_mn100 parax_ee
for i=1:rep
    x_mean(i,:)=mean(parax{i});
    x_var(i,:)=var(parax{i});
end

xi_test=parax{1};
xi_sim=XI{1};
xi_ture=xi_sim{1};

