%function to estimate xi_y for all data files
% consider remove five cosest days.

function [paray,flagy]=predict_logisticy_lag_fun(data,m,grid_length,cov_y)
% cov1_x=cov_x(:,:,i);
% cov1_y=cov_y(:,:,i);
T=1; %number of grid.
grid=T/grid_length/2:T/grid_length:T;

%     [V_x1,D_x1]=eig(cov_x/grid_length);
[V_y1,D_y1]=eig(cov_y/grid_length);

%     D_x=diag(D_x1);
D_y=diag(D_y1);
%     V_x=grid_length^0.5*V_x1;
V_y=grid_length^0.5*V_y1;
V_y(:,end-1)=(2*double(V_y(13,end-1)>0)-1)*V_y(:,end-1);
V_y(:,end)=(2*double(V_y(1,end)>0)-1)*V_y(:,end);


%% get the closest grid of data pair of one subject one day
[ndata,~]=size(data);
for i=1:ndata  % find the closest grid of time point
    [~,grid_id(i)]=min(abs(data(i,3)-grid));
end
data_expand=[data,grid_id'];



%% predict y

parfor j=1:m
    day_close_ind=data_expand(data_expand(:,2)==j,4);
    xi_yd1(j)=sum(V_y(day_close_ind,end)); % first order derivative of xi^z_k,1*m
    xi_yd2(j)=sum(V_y(day_close_ind,end-1));
    xi_yd3(j)=sum(V_y(day_close_ind,end-2));
end
xi_yd=[xi_yd1; xi_yd2; xi_yd3 ];%3*m

covy_diag=diag(cov_y);
% covy_grid_all=exp(covy_diag(grid_id)/2); %N*1

parfor j=1:m
%     sigma_grid=covy_diag(data_expand(data_expand(:,2)~=j,4))/2;
    grid_id1=data_expand((data_expand(:,2)<j-5) | (data_expand(:,2)>j+5)|(data_expand(:,2)==j),4);
    covy_grid1=covy_diag(grid_id1)/2; %grid in denominator
    
%     grid_id2=data_expand((data_expand(:,2)<j-5) | (data_expand(:,2)>j+5),4);
%     covy_grid2=covy_diag(grid_id2)/2; 
    
    phiy_grid=[V_y(grid_id1,end) V_y(grid_id1,end-1) V_y(grid_id1,end-2)]; %estimated phi_x on time grid. dim: kpx*N
    mj=sum(([1:m]<j-5) | ([1:m]>j+5));
    %         xi_y0=[random('Normal',0,D_y(end).^0.5 ,m,1);random('Normal',0,D_y(end-1).^0.5 ,m,1);random('Normal',0,D_y(end-2).^0.5 ,m,1)]; %3m*1,forgot why use this
    xi_y0=[random('Normal',0,D_y(end).^0.5 ,1,1), random('Normal',0,D_y(end-1).^0.5 ,1,1), random('Normal',0,D_y(end-2).^0.5 ,1,1)]; %3*1
    options=optimset('Display','off');
    [paray(j,:),fvaly,flagy(j),outputy]=fminsearch(@logisticy_fun,xi_y0,options,mj,phiy_grid,xi_yd(:,j),covy_grid1);
end

% case with rid the adjacent days
function [logLy]=logisticy_fun(xi_y,mj,phiy_grid,xi_yd,covy_grid1)
% xi_x=xi_x0;
% L1=xi_y*xi_yd +sum(covy_grid2)+log(mj)*length(covy_grid2);
L1=xi_y*xi_yd;
L2=log(exp(xi_y*phiy_grid')+mj*exp(covy_grid1)');%1*N
logLy=-L1+sum(L2);

% case without rid the adjacent days
% function [logLy]=logisticy_fun(xi_y,m,phiy_grid,xi_yd,covy_grid,sigma_grid)
% % xi_x=xi_x0;
% L1=xi_y*xi_yd +sum(sigma_grid);
% L2=log(exp(xi_y*phiy_grid')+m*covy_grid');%1*N
% logLy=-L1+sum(L2);

%% x case
% parfor i=1:n
%     sigma_grid=covx_diag(data_expand(data_expand(:,1)~=i,4))/2;
%     xi_x0=[random('Normal',0,D_x(end).^0.5 ,1,1), random('Normal',0,D_x(end-1).^0.5 ,1,1),  random('Normal',0,D_x(end-2).^0.5 ,1,1)]; %2m*1
%     options=optimset('Display','off');
%     [parax(i,:),fvalx,flagx(i),outputx]=fminsearch(@logisticx_fun,xi_x0,options,n,phix_grid,xi_xd(:,i),covx_grid,sigma_grid);
% end
% function [logLx]=logisticx_fun(xi_x,n,phix_grid,xi_xd,covx_grid,sigma_grid)
% % xi_x=xi_x0;
%
% L1=xi_x*xi_xd +sum(sigma_grid);
% L2=log(exp(xi_x*phix_grid')+(n-1)*covx_grid');%1*N
% logLx=-L1+sum(L2);












