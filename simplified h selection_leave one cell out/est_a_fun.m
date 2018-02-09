%randb_sample:first column is subject when est b;first two columns are subjects
%when est c


function [A_leave,A]=est_a_fun(data,m,n,h,grid_length,randa_sample)

% clear all
% load data100_n100_m100.mat
% data=DATA{1};
% h=0.01;grid_length=50;m=100;n=m;
% na_sample=100;
% randa_sample=sortrows(floor(rand(na_sample,2)*100+1));
% nb_sample=10000;
% randb_sample=sortrows(floor(rand(nb_sample,3)*100+1));

T=1;
% tic
grid=[T/grid_length/2:T/grid_length:T];
na_sample=length(randa_sample(:,1));
% nb_sample=length(randb_sample(:,1));

%get general delta delta1
% tic
parfor i=1:n
%     i=1
    sub_data=data(data(:,1)==i,:);
   for j=1:m
       cell_data=sub_data(sub_data(:,2)==j,3); % get data in each cell
       delta(i,j,:)=sum(exp(-0.5*(bsxfun(@minus, cell_data',grid')/h).^2)/sqrt(2*pi)/h,2);%get delta_ij in each cell
       delta1_tmp=exp(-0.5*(bsxfun(@minus, cell_data',grid')/h).^2)/sqrt(2*pi)/h;%gird_length*#N_ij
       delta1(:,:,i,j)=delta1_tmp*delta1_tmp'; %T2 in A term
   end
end
% toc

%get full term of t1 and t3
% delta_sumj=squeeze(sum(delta,2));%n*ngrid
% tic
parfor i=1:n
    delta_t1_i(:,:,i)=squeeze(delta(i,:,:))' * squeeze(delta(i,:,:));%T1 of ith sub
%    delta_t3_i(:,:,i)=bsxfun(@times, delta_sumj(i,:),delta_sumj(i,:)');
end
% toc
t1a1=sum(delta_t1_i,3); %gird_length*grid_lengh


%leave k,g out for estimation A
t2a1=sum(sum(delta1,3),4); %quantities in T2 of A
t2a2=sum(delta1,3);
t2a3=sum(delta1,4);
% tic
parfor k=1:na_sample
%     i=1
    t1a2=squeeze(delta(:,randa_sample(k,2),:));%n*grid_length
    t1a3=squeeze(delta(randa_sample(k,1),:,:));%grid_length*m
    t1a4=squeeze(delta(randa_sample(k,1),randa_sample(k,2),:)); %grid_length*1
    t1a_k(:,:,k)=t1a1-t1a2'*t1a2-t1a3'*t1a3+t1a4*t1a4';
    t2a_k(:,:,k)=t2a1-t2a2(:,:,:,randa_sample(k,2))-t2a3(:,:,randa_sample(k,1),:)+delta1(:,:,randa_sample(k,1),randa_sample(k,2));
end
% toc

A_leave=t1a_k-t2a_k;
% B_leave=t3b_k-t1b_k;
% C_leave=t4c_k-t1c_k;
A=t1a1-t2a1;
% B=t3b1-t1a1;
% C=t4c1-t1a1;
% toc




