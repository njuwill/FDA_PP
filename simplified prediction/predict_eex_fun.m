function [parax,flagx]=predict_eex_fun(data,n,grid_length,cov_x)

T=1; %number of grid.
grid=T/grid_length/2:T/grid_length:T;
covx_diag=diag(cov_x);
phix_grid_sum=sum(phix_grid_all,1);
[ndata,~]=size(data);
for i=1:ndata  % find the closest grid of time point
    [~,grid_id(i)]=min(abs(data(i,3)-grid));
end
data_expand=[data,grid_id'];

parfor i=1:n
    sub_grid=data_expand(data_expand(:,1)==i,4);
    phix_grid=[V_x(sub_grid,end) V_x(sub_grid,end-1) V_x(sub_grid,end-2)]; %estimated phi_x on time grid. dim: kpx*N    
    covx_grid=exp(covx_diag(sub_grid)/2); %N*1   
    xi_x0=[random('Normal',0,D_x(end).^0.5 ,1,1), random('Normal',0,D_x(end-1).^0.5 ,1,1),  random('Normal',0,D_x(end-2).^0.5 ,1,1)]; %2m*1
    options=optimset('Display','off');
    [parax(i,:),fvalx,flagx(i),outputx]=fsolve(@eex_fun,xi_x0,options,n,phix_grid,xi_xd(:,i),covx_grid,sigma_grid);
end
function [Ux]=eex_fun(xi_x,n_all,phix_grid,covx_grid,phix_grid_sum)
% xi_x=xi_x0;

U1=phix_grid.*kron(ones(1,3),exp(covx_grid/2)./exp(phix_grid*xi_x'))*(n_all-1);
Ux=sum(U1,1)-phix_grid_sum+sum(phix_grid,1);


