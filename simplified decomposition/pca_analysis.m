
clear all

grid_length=100; 
T=1;grid=[T/grid_length/2:T/grid_length:T];
% h=0.01;
rep=300;

load cov_h002_sim300_mn400
cov_x=cov_xb;cov_y=cov_yb;cov_z=cov_zb;
% cov_x=cov_xs;cov_y=cov_ys;cov_z=cov_zs;

for i=1:rep
     [V_x,D_x]=eig(cov_x(:,:,i)/grid_length);  
     [V_y,D_y]=eig(cov_y(:,:,i)/grid_length); 
     [V_z,D_z]=eig(cov_z(:,:,i)/grid_length);   
     Dx_all(:,i)=diag(D_x);   Dy_all(:,i)=diag(D_y);  Dz_all(:,i)=diag(D_z);        
     Vx_all(:,:,i)=grid_length^0.5*V_x;  
     Vy_all(:,:,i)=grid_length^0.5*V_y;  
     Vz_all(:,:,i)=grid_length^0.5*V_z;
     [Vx_b(:,:,i), Vy_b(:,:,i),Vz_b(:,:,i)]=change_dir(grid_length,Vx_all(:,:,i),Vy_all(:,:,i),Vz_all(:,:,i));
end
save DVb_mn400 Dx_all Dy_all Dz_all Vx_b Vy_b Vz_b
save DVs_mn400 Dx_all Dy_all Dz_all Vx_s Vy_s Vz_s


%plot boxplot of eigenvalues
% figure
% subplot(3,2,1);boxplot(squeeze(Dx_all(end,:,:)))
% subplot(3,2,3);boxplot(squeeze(Dy_all(end,:,:)))
% subplot(3,2,4);boxplot(squeeze(Dy_all(end-1,:,:)))
% subplot(3,2,5);boxplot(squeeze(Dz_all(end,:,:)))
% subplot(3,2,6);boxplot(squeeze(Dz_all(end-1,:,:)))
%get std of eigenvalues

load DVb_mn400
mean_all=[mean(squeeze(Dx_all(end,:,:)))',mean(squeeze(Dx_all(end-1,:,:)))',...
mean(squeeze(Dy_all(end,:,:)))',mean(squeeze(Dy_all(end-1,:,:)))'...
mean(squeeze(Dz_all(end,:,:)))',mean(squeeze(Dz_all(end-1,:,:)))']';

std_all=[std(squeeze(Dx_all(end,:,:)))',std(squeeze(Dx_all(end-1,:,:)))',...
 std(squeeze(Dy_all(end,:,:)))',std(squeeze(Dy_all(end-1,:,:)))',...
 std(squeeze(Dz_all(end,:,:)))' ,std(squeeze(Dz_all(end-1,:,:)))']';
[mean_all std_all]

Vx_all=Vx_b;Vy_all=Vy_b;Vz_all=Vz_b;
% Vx_all=Vx_s;Vy_all=Vy_s;Vz_all=Vz_s;
figure
subplot(2,3,1);plot(grid,squeeze(Vx_all(:,end,:)));ylim([-2,2])
subplot(2,3,4);plot(grid,squeeze(Vx_all(:,end-1,:)));ylim([-2,2])

subplot(2,3,2);plot(grid,squeeze(Vy_all(:,end,:)));ylim([-2,2])
subplot(2,3,5);plot(grid,squeeze(Vy_all(:,end-1,:)));ylim([-2,2])

subplot(2,3,3);plot(grid,squeeze(Vz_all(:,end,:)));ylim([-2,2])
subplot(2,3,6);plot(grid,squeeze(Vz_all(:,end-1,:)));ylim([-2,2])


%plot eigenfunction of one sample size
Vx_all=Vx_b;Vy_all=Vy_b;Vz_all=Vz_b;
% Vx_all=Vx_s;Vy_all=Vy_s;Vz_all=Vz_s;
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

% plot eigenfunction of two sample size
load DVb_mn200 
Vx_all200=Vx_b;Vy_all200=Vy_b;Vz_all200=Vz_b;
load DVb_mn400 
Vx_all400=Vx_b;Vy_all400=Vy_b;Vz_all400=Vz_b;
subplot(3,4,1);plot(squeeze(Vx_all200(:,end,:)));ylim([-2,2])
subplot(3,4,2);plot(squeeze(Vx_all200(:,end-1,:)));ylim([-2,2])
subplot(3,4,3);plot(squeeze(Vx_all400(:,end,:)));ylim([-2,2])
subplot(3,4,4);plot(squeeze(Vx_all400(:,end-1,:)));ylim([-2,2])

subplot(3,4,5);plot(squeeze(Vy_all200(:,end,:)));ylim([-2,2])
subplot(3,4,6);plot(squeeze(Vy_all200(:,end-1,:)));ylim([-2,2])
subplot(3,4,7);plot(squeeze(Vy_all400(:,end,:)));ylim([-2,2])
subplot(3,4,8);plot(squeeze(Vy_all400(:,end-1,:)));ylim([-2,2])

subplot(3,4,9);plot(squeeze(Vz_all200(:,end,:)));ylim([-2,2])
subplot(3,4,10);plot(squeeze(Vz_all200(:,end-1,:)));ylim([-2,2])
subplot(3,4,11);plot(squeeze(Vz_all400(:,end,:)));ylim([-2,2])
subplot(3,4,12);plot(squeeze(Vz_all400(:,end-1,:)));ylim([-2,2])
