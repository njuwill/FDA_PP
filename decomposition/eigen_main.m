clear all
load('C:\Users\mwang\Dropbox\finance project\pegasus\data\sim300\data300_nm100_rho5.mat')
% load data100_nm200_rho5.mat
grid_length=100;
T=1;
rep=300;
m=100;n=m;
tic
load h100
for i=1:rep
    i
    data=DATA{i};
    %     hvec=[ha_200(i) hb_200(i) hc_200(i) hrho_200(i)];
    hvec=[ha_100(i) hb_100(i) hc_100(i) hrho_100(i)];
    [cov_x(:,:,i),cov_y(:,:,i),cov_z(:,:,i)]=get_cov_fun(data,m,n,hvec,grid_length);
end
toc
save cov_sepH_sim300_mn100 cov_x cov_y cov_z

%case nm200
load('C:\Users\mwang\Dropbox\finance project\pegasus\data\sim300\data300_nm200_rho5.mat')
% load data100_nm200_rho5.mat
grid_length=100;
T=1;
rep=300;
m=200;n=m;
tic
load h200
for i=1:rep
    i
    data=DATA{i};
    hvec=[ha_200(i) hb_200(i) hc_200(i) hrho_200(i)];
%     hvec=[ha_100(i) hb_100(i) hc_100(i) hrho_100(i)];
    [cov_x(:,:,i),cov_y(:,:,i),cov_z(:,:,i)]=get_cov_fun(data,m,n,hvec,grid_length);
end
toc
save cov_sepH_sim300_mn200 cov_x cov_y cov_z


%% summary
clear all
load cov_sepH_sim300_mn200
for i=1:rep
    [V_x,D_x]=eig(cov_x(:,:,i)/grid_length);
    [V_y,D_y]=eig(cov_y(:,:,i)/grid_length);
    [V_z,D_z]=eig(cov_z(:,:,i)/grid_length);
    Dx_all(:,i)=diag(D_x);   Dy_all(:,i)=diag(D_y);  Dz_all(:,i)=diag(D_z);
    Vx_all(:,:,i)=grid_length^0.5*V_x;
    Vy_all(:,:,i)=grid_length^0.5*V_y;
    Vz_all(:,:,i)=grid_length^0.5*V_z;
end
[mean(squeeze(Dx_all(end,:,:)))',std(squeeze(Dx_all(end,:,:)))',...
    mean(squeeze(Dy_all(end,:,:)))',std(squeeze(Dy_all(end,:,:)))',...
    mean(squeeze(Dy_all(end-1,:,:)))',std(squeeze(Dy_all(end-1,:,:)))'...
    mean(squeeze(Dz_all(end,:,:)))',std(squeeze(Dz_all(end,:,:)))' ...
    mean(squeeze(Dz_all(end-1,:,:)))',std(squeeze(Dz_all(end-1,:,:)))'...
    mean(squeeze(Dz_all(end-2,:,:)))',std(squeeze(Dz_all(end-2,:,:)))']'


figure
subplot(3,3,1);plot(squeeze(Vx_all(:,end,:)));ylim([-2,2])
subplot(3,3,2);plot(squeeze(Vx_all(:,end-1,:)));ylim([-2,2])
subplot(3,3,3);plot(squeeze(Vx_all(:,end-2,:)));ylim([-2,2])
subplot(3,3,4);plot(squeeze(Vy_all(:,end,:)));ylim([-2,2])
subplot(3,3,5);plot(squeeze(Vy_all(:,end-1,:)));ylim([-2,2])
% Vy2=kron(squeeze((2*double(Vy_all(5,end-1,:)>0)-1))',ones(20,1)).*squeeze(Vy_all(:,end-1,:));
% subplot(3,3,5);plot(squeeze(Vy2));ylim([-2,2])
subplot(3,3,6);plot(squeeze(Vy_all(:,end-2,:)));ylim([-2,2])
% Vz1=kron(squeeze((2*double(Vz_all(3,end,:)>0)-1))',ones(20,1)).*squeeze(Vz_all(:,end,:));
% subplot(3,3,7);plot(squeeze(Vz1));ylim([-2,2])
% Vz2=kron(squeeze((2*double(Vz_all(1,end-1,:)>0)-1))',ones(20,1)).*squeeze(Vz_all(:,end-1,:));
% subplot(3,3,8);plot(squeeze(Vz2));ylim([-2,2])
subplot(3,3,7);plot(squeeze(Vz_all(:,end,:)));ylim([-2,2])
subplot(3,3,8);plot(squeeze(Vz_all(:,end-1,:)));ylim([-2,2])
subplot(3,3,9);plot(squeeze(Vz_all(:,end-2,:)));ylim([-2,2])
