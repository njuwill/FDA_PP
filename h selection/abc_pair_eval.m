%get the leave one out pair when estimate X^(-k)(u,v;h)
%then find the cloest grid of each pair, renturn the grid index

function [paira,pairb,pairc,paira_same]=abc_pair_eval(data,m,n,grid_length,randa_sample,randb_sample)
% %test input
% load data100_n100_m100
% data=DATA{1};
% N=100;M=N;
% grid_length=50;
% na_sample=100;
% randa_sample=sortrows(floor(rand(na_sample,2)*100+1));
% nb_sample=1000000;
% randb_sample=sortrows(floor(rand(nb_sample,3)*100+1));

na_sample=length(randa_sample(:,1));
nb_sample=length(randb_sample(:,1));
[ndata,~]=size(data);
T=1;
grid=T/grid_length/2:T/grid_length:T;
for i=1:ndata  % find the closest grid of time point
    %grid_id(i)=find(abs(data(i,3)-grid)==min(abs(data(i,3)-grid)),1); %if two exist,take min
    [~,grid_id(i)]=min(abs(data(i,3)-grid));
end
data=[data,grid_id'];

parfor i=1:n
    sub_data=data(data(:,1)==i,:);
    for j=1:m
        cell_data{i,j}=sub_data(sub_data(:,2)==j,4); % get data in each cell
    end
end

% get the cloest grid for a: get all first, then get rid of pair of same point
parfor k=1:na_sample
        [pair_cell1,pair_cell2]=meshgrid(cell_data{randa_sample(k,1),randa_sample(k,2)},cell_data{randa_sample(k,1),randa_sample(k,2)});
        [pair_a1,pair_a2]=deal(reshape(pair_cell1,[],1),reshape(pair_cell2,[],1));
        paira{k}=sub2ind([grid_length grid_length],pair_a1,pair_a2);
        pair_cell_same=[cell_data{randa_sample(k,1),randa_sample(k,2)},cell_data{randa_sample(k,1),randa_sample(k,2)}];
        paira_same{k}=sub2ind([grid_length grid_length],pair_cell_same(:,1),pair_cell_same(:,2));
end

parfor k=1:nb_sample
    [pair_sub1,pair_sub2]=meshgrid(cell_data{randb_sample(k,1),randb_sample(k,2)},cell_data{randb_sample(k,1),randb_sample(k,3)}); 
    [pair_b1,pair_b2]=deal(reshape(pair_sub1,[],1),reshape(pair_sub2,[],1));
    pairb{k}= sub2ind([grid_length grid_length],pair_b1,pair_b2);
end


% tic
parfor k=1:nb_sample
    [pair_day1,pair_day2]=meshgrid(cell_data{randb_sample(k,1),randb_sample(k,3)},cell_data{randb_sample(k,2),randb_sample(k,3)}); 
    [pair_c1,pair_c2]=deal(reshape(pair_day1,[],1),reshape(pair_day2,[],1));
    pairc{k}= sub2ind([grid_length grid_length],pair_c1,pair_c2);
end
% toc
