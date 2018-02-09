%get the data pairs of same person same day  when estimate lambda_sb(u,v)
%then find the cloest grid of each pair, renturn the grid index

%psi_zb: grid_length*2
% function [pair_sumj,pair_ij,pair_count]=pairA_idx(data_b,data_s,n,m,grid_length)
function [pair_ij,f2]=pair_eez_fun(data_b,data_s,V_zb,V_zs,n_covz,n,m,grid_length)

% % % %test input
% load data100_mn100_bsxycorz
% data_b=DATA_B{1};
% data_s=DATA_S{1};
% n=100;m=n;
% grid_length=50;

psi_zb=fliplr(V_zb(:,end-n_covz+1:end)); %50 *2 , first two eigenfunctions of buying proceess
psi_zs=fliplr(V_zs(:,end-n_covz+1:end)); %50 *2, first two eigenfunctions of selling proceess


[ndata_b,~]=size(data_b);
[ndata_s,~]=size(data_s);
T=1;
grid=T/grid_length/2:T/grid_length:T;
clear grid_idb grid_ids
for i=1:ndata_b  % find the closest grid of time point
    %grid_id(i)=find(abs(data(i,3)-grid)==min(abs(data(i,3)-grid)),1); %if two exist,take min
    [~,grid_idb(i)]=min(abs(data_b(i,3)-grid));
end
data_b=[data_b,grid_idb'];
for i=1:ndata_s  % find the closest grid of time point
    %grid_id(i)=find(abs(data(i,3)-grid)==min(abs(data(i,3)-grid)),1); %if two exist,take min
    [~,grid_ids(i)]=min(abs(data_s(i,3)-grid));
end
data_s=[data_s,grid_ids'];


% get sum over phi_z on diff persons, diff days
for i=1:n
%          i
    sub_data_b=data_b(data_b(:,1)==i,:);
    sub_data_s=data_s(data_s(:,1)==i,:);
    for j=1:m
        cell_data_b=sub_data_b(sub_data_b(:,2)==j,4); % get data in each cell
        cell_data_s=sub_data_s(sub_data_s(:,2)==j,4);
        cell_phi_b(i,j,:)=sum(psi_zb(cell_data_b,:),1);
        cell_phi_s(i,j,:)=sum(psi_zs(cell_data_s,:),1);
        cell_multiply_tmp(:,:,i,j)=sum(psi_zb(cell_data_b,:),1)'*sum(psi_zs(cell_data_s,:),1);
    end
end
% 
% for i=1:n
%     for j=1:m
%         cell_multiply_tmp(:,:,i,j)=squeeze(cell_phi_b(i,j,:))*squeeze(cell_phi_s(i,j,:))';
%     end
% end

f2=squeeze(sum(sum(cell_phi_b,1)))*squeeze(sum(sum(cell_phi_s)))'-...
    squeeze(sum(cell_phi_b))'*squeeze(sum(cell_phi_s))-...
    squeeze(sum(cell_phi_b,2))'*squeeze(sum(cell_phi_s,2))+...
    squeeze(sum(sum(cell_multiply_tmp,3),4));




%% get pair (t_b,t_s) from same person, same day(7s)
% tic
pair_ij=[];
for i=1:n
    sub_data_b=data_b(data_b(:,1)==i,:);
    sub_data_s=data_s(data_s(:,1)==i,:);
%     [pair_cell1,pair_cell2]=meshgrid(sub_data_b(:,4),sub_data_s(:,4)');
%     [pair_a1,pair_a2]=deal(reshape(pair_cell1,[],1),reshape(pair_cell2,[],1));
%     pair_tmp=sub2ind([grid_length grid_length],pair_a1,pair_a2);
%     pair_index_sumj(i,:)=sum(bsxfun(@eq,pair_tmp,(1:2500))); %1*2500  
    for j=1:m
        % get the cloest grid for a: get all first, then get pair of two processes
        cell_data_b=sub_data_b(sub_data_b(:,2)==j,4); % get data in each cell
        cell_data_s=sub_data_s(sub_data_s(:,2)==j,4); % get data in each cell
        [pair_cell1,pair_cell2]=meshgrid(cell_data_b,cell_data_s');
        [pair_a1,pair_a2]=deal(reshape(pair_cell1,[],1),reshape(pair_cell2,[],1));
        pair_tmp=[pair_a1,pair_a2];
        pair_ij=[pair_ij;pair_tmp];
    end
end
% toc
