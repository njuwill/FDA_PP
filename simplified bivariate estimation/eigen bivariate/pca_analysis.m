
clear all

grid_length=50; 
T=1;grid=[T/grid_length/2:T/grid_length:T];
rep=300;

load cov_h002_sim300_mn400
% estimate eigenvalues and eigenfuctions without lag info.
cov_x=cov_xb;cov_y=cov_yb;cov_z=cov_zb;
cov_x=cov_xs;cov_y=cov_ys;cov_z=cov_zs;
for i=1:rep
     [V_x,D_x]=eig(cov_x(:,:,i)/grid_length);  
     [V_y,D_y]=eig(cov_y(:,:,i)/grid_length); 
     [V_z,D_z]=eig(cov_z(:,:,i)/grid_length);   
     Dx_all(:,i)=diag(D_x);   Dy_all(:,i)=diag(D_y);  Dz_all(:,i)=diag(D_z);        
     Vx_all(:,:,i)=grid_length^0.5*V_x;  
     Vy_all(:,:,i)=grid_length^0.5*V_y;  
     Vz_all(:,:,i)=grid_length^0.5*V_z;
     
end



%get mean and std of eigenvalues
mean_all=[mean(squeeze(Dx_all(end,:,:)))',mean(squeeze(Dx_all(end-1,:,:)))',...
mean(squeeze(Dy_all(end,:,:)))',mean(squeeze(Dy_all(end-1,:,:)))'...
mean(squeeze(Dz_all(end,:,:)))',mean(squeeze(Dz_all(end-1,:,:)))']';

std_all=[std(squeeze(Dx_all(end,:,:)))',std(squeeze(Dx_all(end-1,:,:)))',...
 std(squeeze(Dy_all(end,:,:)))',std(squeeze(Dy_all(end-1,:,:)))',...
 std(squeeze(Dz_all(end,:,:)))' ,std(squeeze(Dz_all(end-1,:,:)))']';
[mean_all std_all]

%plot eigenfunction 
figure
Vx1=kron(squeeze((2*double(Vx_all(50,end,:)>0)-1))',ones(100,1)).*squeeze(Vx_all(:,end,:));
subplot(3,3,1);plot(squeeze(Vx1));ylim([-2,2])
% subplot(3,3,1);plot(squeeze(Vx_all(:,end,:)));ylim([-2,2])
Vx2=kron(squeeze((2*double(Vx_all(1,end-1,:)>0)-1))',ones(100,1)).*squeeze(Vx_all(:,end-1,:));
subplot(3,3,2);plot(squeeze(Vx2));ylim([-2,2])
% subplot(3,3,2);plot(squeeze(Vx_all(:,end-1,:)));ylim([-2,2])
subplot(3,3,3);plot(squeeze(Vx_all(:,end-2,:)));ylim([-2,2])
subplot(3,3,4);plot(squeeze(Vy_all(:,end,:)));ylim([-2,2])
% subplot(3,3,5);plot(squeeze(Vy_all(:,end-1,:)));ylim([-2,2])
Vy2=kron(squeeze((2*double(Vy_all(25,end-1,:)>0)-1))',ones(100,1)).*squeeze(Vy_all(:,end-1,:));
subplot(3,3,5);plot(squeeze(Vy2));ylim([-2,2])
subplot(3,3,6);plot(squeeze(Vy_all(:,end-2,:)));ylim([-2,2])
Vz1=kron(squeeze((2*double(Vz_all(50,end,:)>0)-1))',ones(100,1)).*squeeze(Vz_all(:,end,:));
subplot(3,3,7);plot(squeeze(Vz1));ylim([-2,2])
Vz2=kron(squeeze((2*double(Vz_all(13,end-1,:)>0)-1))',ones(100,1)).*squeeze(Vz_all(:,end-1,:));
subplot(3,3,8);plot(squeeze(Vz2));ylim([-2,2])
% subplot(3,3,7);plot(squeeze(Vz_all(:,end,:)));ylim([-2,2])
% subplot(3,3,8);plot(squeeze(Vz_all(:,end-1,:)));ylim([-2,2])
subplot(3,3,9);plot(squeeze(Vz_all(:,end-2,:)));ylim([-2,2])


%plot boxplot of eigenvalues
% figure
% subplot(3,2,1);boxplot(squeeze(Dx_all(end,:,:)))
% subplot(3,2,2);boxplot(squeeze(Dx_all(end-1,:,:)))
% subplot(3,2,3);boxplot(squeeze(Dy_all(end,:,:)))
% subplot(3,2,4);boxplot(squeeze(Dy_all(end-1,:,:)))
% subplot(3,2,5);boxplot(squeeze(Dz_all(end,:,:)))
% subplot(3,2,6);boxplot(squeeze(Dz_all(end-1,:,:)))


