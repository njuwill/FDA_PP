clear all

load('D:\Ming Wang\sim data\sim300final\data300_mn400_bs.mat')

grid_length=100;
T=1;
rep=300;
m=400;n=m;
tic
% load h100
for i=1:300
    i
    datab=DATA_B{i};
    datas=DATA_S{i};
    %     hvec=[ha_200(i) hb_200(i) hc_200(i) hrho_200(i)];
%     hvec=[ha_100(i) hb_100(i) hc_100(i) hrho_100(i)];
    h=0.02;
    [q_x(:,:,i),q_y(:,:,i),q_z(:,:,i)]=get_covbs_fun(datab,datas,m,n,h,grid_length);
end
toc
save q_h002_mn400 q_x q_y q_z

% load DVb_mn200
% load DVs_mn200

load DVb_mn400
load DVs_mn400

T=1;
grid=T/grid_length/2:T/grid_length:T;


clear sigma_x sigma_y sigma_z
parfor i=1:300
    sigma_x(:,:,i)=flip(Vx_b(:,end-1:end,i),2)'*q_x(:,:,i)*flip(Vx_s(:,end-1:end,i),2)*1/grid_length^2;
    sigma_y(:,:,i)=flip(Vy_b(:,end-1:end,i),2)'*q_y(:,:,i)*flip(Vy_s(:,end-1:end,i),2)*1/grid_length^2;
    sigma_z(:,:,i)=flip(Vz_b(:,end-1:end,i),2)'*q_z(:,:,i)*flip(Vz_s(:,end-1:end,i),2)*1/grid_length^2;
end
mean(sigma_x,3)
[std(sigma_x(1,1,:)) std(sigma_x(1,2,:));std(sigma_x(2,1,:)) std(sigma_x(2,2,:))]

mean(sigma_y,3)
[std(sigma_y(1,1,:)) std(sigma_y(1,2,:));std(sigma_y(2,1,:)) std(sigma_y(2,2,:))]

mean(sigma_z,3)
[std(sigma_z(1,1,:)) std(sigma_z(1,2,:));std(sigma_z(2,1,:)) std(sigma_z(2,2,:))]



for i=1:rep
    i=1
     [V_x,D_x]=eig(q_x(:,:,i)/grid_length);  
     [V_y,D_y]=eig(q_y(:,:,i)/grid_length); 
     [V_z,D_z]=eig(q_z(:,:,i)/grid_length);   
     Dx_all(:,i)=diag(D_x);   Dy_all(:,i)=diag(D_y);  Dz_all(:,i)=diag(D_z);        
     Vx_all(:,:,i)=grid_length^0.5*V_x;  
     Vy_all(:,:,i)=grid_length^0.5*V_y;  
     Vz_all(:,:,i)=grid_length^0.5*V_z;
     
end





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
%% test
i=1;
data=DATA{i};
tic
[estA,estB,estC,estD]=estABCD_fun(data,m,n,h,grid_length);
toc
grid_length=100;
[estA1]=estA_fun(data,m,n,h,grid_length);
[estB1]=estB_fun(data,m,n,h,grid_length);
[estC1]=estC_fun(data,m,n,h,grid_length);
[estD1]=estD_fun(data,m,n,h,grid_length);
head((estD-estD1)./estD)
head((estA-estA1))