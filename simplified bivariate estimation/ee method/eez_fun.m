%estimating equation method to estimate gamma of cov(xi_zb,xi_zs)
function [f]=eez_fun(sigma0,n,m,n_covz,f2,pair_ij,cov_x,cov_y,...
    V_xb,V_yb,V_zb,V_xs,V_ys,V_zs)

psi_xb=[V_xb(:,end) V_xb(:,end-1)]; %50 *1 , first two eigenfunctions of buying proceess
psi_xs=[V_xs(:,end) V_xs(:,end-1)]; %50 *1, first two eigenfunctions of selling proceess
psi_yb=[V_yb(:,end) V_yb(:,end-1)]; %50 *2 , first two eigenfunctions of buying proceess
psi_ys=[V_ys(:,end) V_ys(:,end-1)]; %50 *2, first two eigenfunctions of selling proceess
psi_zb=fliplr(V_zb(:,end-n_covz+1:end)); %50 *2 , first two eigenfunctions of buying proceess
psi_zs=fliplr(V_zs(:,end-n_covz+1:end)); %50 *2, first two eigenfunctions of selling proceess


sigma_xbs=cov_x(1:2,3:4);
sigma_ybs=cov_y(1:2,3:4);
sigma_zbs=reshape(sigma0,n_covz,n_covz);


rx_tmp=sum(psi_xb(pair_ij(:,1),:)*sigma_xbs.*psi_xs(pair_ij(:,2),:),2);
ry_tmp=sum(psi_yb(pair_ij(:,1),:)*sigma_ybs.*psi_ys(pair_ij(:,2),:),2);
rz_tmp=sum(psi_zb(pair_ij(:,1),:)*sigma_zbs.*psi_zs(pair_ij(:,2),:),2);
f1=(psi_zb(pair_ij(:,1),:)./(exp(rx_tmp+ry_tmp+rz_tmp)*ones(1,n_covz)))'*psi_zs(pair_ij(:,2),:);

f=(n-1)*(m-1)*f1-f2;

% tic
% n_pair=length(pair_ij(:,1));
% parfor i=1:n_pair
% %     i
%     rx_tmp=exp(psi_xb(pair_ij(i,1),:)*sigma_xbs*psi_xs(pair_ij(i,2),:)');
%     ry_tmp=exp(psi_yb(pair_ij(i,1),:)*sigma_ybs*psi_ys(pair_ij(i,2),:)');
%     f1_array(:,:,i)=psi_zb(pair_ij(i,1),:)'*psi_zs(pair_ij(i,2),:)/rx_tmp/ry_tmp;
% end
% parfor i=1:n_pair
%     rz_array(:,:,i)=ones(n_covz)*exp(psi_zb(pair_ij(i,1),:)*sigma_zbs*psi_zs(pair_ij(i,2),:)');
%     f1(:,:,i)=f1_array(:,:,i)./rz_array(:,:,i);
% end
%  f=(n-1)*(m-1)*sum(f1,3)-f2;
% toc

end
