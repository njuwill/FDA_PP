load data100_n100_m100.mat
data=DATA{2};
hvec=0.0025:0.0025:0.06;
grid_length=100;m=100;n=m;

[pair_cell1,pair_cell2]=meshgrid([1:n],[1:n]');
[pair_a1,pair_a2]=deal(reshape(pair_cell1,[],1),reshape(pair_cell2,[],1));
randa_sample=[pair_a1 pair_a2];
% nb_sample=10000*percentage;
nb_sample=10000;
randb_sample=sortrows(floor(rand(nb_sample,3)*100+1));
na_sample=n*n;

T=1;
grid=[T/grid_length/2:T/grid_length:T];

%% get the pair
for i=1:n
    sub_data=data(data(:,1)==i,:);
    for j=1:m
        cell_data{i,j}=sub_data(sub_data(:,2)==j,3); % get data in each cell
    end
end

pair_aall=[];
pair_a_same=[];
parfor k=1:na_sample
    [pair_cell1,pair_cell2]=meshgrid(cell_data{randa_sample(k,1),randa_sample(k,2)},cell_data{randa_sample(k,1),randa_sample(k,2)});
    [pair_a1,pair_a2]=deal(reshape(pair_cell1,[],1),reshape(pair_cell2,[],1));
    pair_atmp=[pair_a1,pair_a2];
    pair_aall=[pair_aall;pair_atmp];
    pair_cell_same=[cell_data{randa_sample(k,1),randa_sample(k,2)},cell_data{randa_sample(k,1),randa_sample(k,2)}];    
    pair_a_same=[pair_a_same;pair_cell_same];
end

pair_ball=[];
parfor k=1:nb_sample
    [pair_sub1,pair_sub2]=meshgrid(cell_data{randb_sample(k,1),randb_sample(k,2)},cell_data{randb_sample(k,1),randb_sample(k,3)}); 
    [pair_b1,pair_b2]=deal(reshape(pair_sub1,[],1),reshape(pair_sub2,[],1));
    pair_btmp=[pair_b1,pair_b2];
    pair_ball=[pair_ball;pair_btmp]
end


pair_call=[];
parfor k=1:nb_sample
    [pair_day1,pair_day2]=meshgrid(cell_data{randb_sample(k,1),randb_sample(k,3)},cell_data{randb_sample(k,2),randb_sample(k,3)}); 
    [pair_c1,pair_c2]=deal(reshape(pair_day1,[],1),reshape(pair_day2,[],1));
    pair_ctmp=[pair_c1,pair_c2];
    pair_call=[pair_call;pair_ctmp]
end
%% count the number
pair_all=pair_call;
h=0.01
parfor i=1:grid_length
    for j=1:grid_length
        num_logic(i,j)=sum(double(abs(pair_all(:,1)-ones(length(pair_all),1)*grid(i))<h).*...
            double(abs(pair_all(:,2)-ones(length(pair_all),1)*grid(j))<h));
    end
end
sum(sum(num_logic))/10000

% for a
pair_all1=pair_aall;
pair_all2=pair_a_same;
h=0.01
parfor i=1:grid_length
    for j=1:grid_length
        num_logic(i,j)=sum(double(abs(pair_all1(:,1)-ones(length(pair_all1),1)*grid(i))<h).*...
            double(abs(pair_all1(:,2)-ones(length(pair_all1),1)*grid(j))<h));
             num_logic1(i,j)=sum(double(abs(pair_all2(:,1)-ones(length(pair_all2),1)*grid(i))<h).*...
            double(abs(pair_all2(:,2)-ones(length(pair_all2),1)*grid(j))<h));
    end
end
sum(sum(num_logic-num_logic1))/10000
