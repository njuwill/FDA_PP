%randb_sample:first column is subject when est b;first two columns are subjects
%when est c


function [A_leave,B_leave,C_leave,A,B,C]=est_abc_fun(data,m,n,h,grid_length,randa_sample,randb_sample)

% clear all
% load data100_n100_m100.mat
% data=DATA{2};
% h=0.01;grid_length=50;m=100;n=m;
% na_sample=100;
% randa_sample=sortrows(floor(rand(na_sample,2)*100+1));
% nb_sample=10000;
% randb_sample=sortrows(floor(rand(nb_sample,3)*100+1));

T=1;
% tic
grid=[T/grid_length/2:T/grid_length:T];
na_sample=length(randa_sample(:,1));
nb_sample=length(randb_sample(:,1));

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
delta_sumj=squeeze(sum(delta,2));%n*ngrid
% tic
parfor i=1:n
    delta_t1_i(:,:,i)=squeeze(delta(i,:,:))' * squeeze(delta(i,:,:));%T1 of ith sub
   delta_t3_i(:,:,i)=bsxfun(@times, delta_sumj(i,:),delta_sumj(i,:)');
end
% toc
t1a1=sum(delta_t1_i,3); %gird_length*grid_lengh
t3b1=sum(delta_t3_i,3);

delta_sumi=squeeze(sum(delta,1));%m*ngrid
parfor j=1:m
    delta_t4_j(:,:,j)=bsxfun(@times,delta_sumi(j,:),delta_sumi(j,:)'); %T4 of jth day
end
t4c1=sum(delta_t4_j,3);

%leav k,g1,g2 out for estimation B
% tic
parfor k=1:nb_sample
    delta_sumg1g2=squeeze(sum(delta(:,randb_sample(k,2:3),:),2));%n*ngrid
    delta_sumg1g2_k=squeeze(sum(delta(randb_sample(k,1),randb_sample(k,2:3),:),2));%ngrid*1
    t3b2=delta_sumj'*delta_sumg1g2;
    t3b3=delta_sumg1g2'*delta_sumg1g2;
    t3b4=delta_sumj(randb_sample(k,1),:)'*delta_sumj(randb_sample(k,1),:);
    t3b5=delta_sumj(randb_sample(k,1),:)'*delta_sumg1g2_k';
    t3b6=delta_sumg1g2_k*delta_sumg1g2_k';
    t3b_k(:,:,k)=t3b1-t3b2-t3b2'+t3b3-t3b4+t3b5+t3b5'-t3b6;
    
    t1b1=t1a1;
    t1b2_g1=squeeze(delta(:,randb_sample(k,2),:));%n*grid_length
    t1b2_g2=squeeze(delta(:,randb_sample(k,3),:));%n*grid_length
    t1b2=t1b2_g1'*t1b2_g1+t1b2_g2'*t1b2_g2;
    t1b3=squeeze(delta(randb_sample(k,1),:,:))'*squeeze(delta(randb_sample(k,1),:,:));
    t1b4_g1=t1b2_g1(randb_sample(k,1),:);%1*grid_length
    t1b4_g2=t1b2_g2(randb_sample(k,1),:);
    t1b4=t1b4_g1'*t1b4_g1+t1b4_g2'*t1b4_g2;
    t1b_k(:,:,k)=t1b1-t1b2-t1b3+t1b4;
end
% toc

%leave k1,k2,g out for estimation C
% tic
parfor k=1:nb_sample
    delta_sumk1k2=squeeze(sum(delta(randb_sample(k,1:2),:,:),1));%m*ngrid
    delta_sumk1k2_k=squeeze(sum(delta(randb_sample(k,1:2),randb_sample(k,3),:),1));%ngrid*1
    t4c2=delta_sumi'*delta_sumk1k2;
    t4c3=delta_sumk1k2'*delta_sumk1k2;
    t4c4=delta_sumi(randb_sample(k,3),:)'*delta_sumi(randb_sample(k,3),:);
    t4c5=delta_sumi(randb_sample(k,3),:)'*delta_sumk1k2_k';
    t4c6=delta_sumk1k2_k*delta_sumk1k2_k';
    t4c_k(:,:,k)=t4c1-t4c2-t4c2'+t4c3-t4c4+t4c5+t4c5'-t4c6;
    
    t1c1=t1a1;
    t1c2_g1=squeeze(delta(randb_sample(k,1),:,:));%m*grid_length
    t1c2_g2=squeeze(delta(randb_sample(k,2),:,:));%m*grid_length
    t1c2=t1c2_g1'*t1c2_g1+t1c2_g2'*t1c2_g2;
    t1c3=squeeze(delta(:,randb_sample(k,3),:))'*squeeze(delta(:,randb_sample(k,3),:));
    t1c4_g1=t1c2_g1(randb_sample(k,3),:);%1*grid_length
    t1c4_g2=t1c2_g2(randb_sample(k,3),:);
    t1c4=t1c4_g1'*t1c4_g1+t1c4_g2'*t1c4_g2;
    t1c_k(:,:,k)=t1c1-t1c2-t1c3+t1c4;
end
% toc



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
B_leave=t3b_k-t1b_k;
C_leave=t4c_k-t1c_k;
A=t1a1-t2a1;
B=t3b1-t1a1;
C=t4c1-t1a1;
% toc
%% rho rho_leave estimation
% parfor i=1:n
% %     i=1
%     sub_data=data(data(:,1)==i,:);
%    for j=1:m
%        cell_data=sub_data(sub_data(:,2)==j,3); % get data in each cell
%        delta1_tmp=exp(-0.5*(bsxfun(@minus, cell_data',grid')/h).^2)/sqrt(2*pi)/h;%gird_length*#N_ij
%        delta2(:,i,j)=sum(delta1_tmp,2);      
%    end
% end
% 
% rho1=sum(sum(delta2,2),3);
% 
% parfor k=1:na_sample
%    rho2=sum(delta2(:,:,randa_sample(k,2)),2);%n*grid_length
%     rho3=sum(delta2(:,randa_sample(k,1),:),3);%n*grid_length
%     rho4=delta2(:,randa_sample(k,1),randa_sample(k,1));%n*grid_length
%     rho_leave(:,k)=rho1-rho2-rho3+rho4;
% end
% rho=rho1;

% rho=sum(exp(-0.5*(bsxfun(@minus, data(:,3),grid)/h).^2)/sqrt(2*pi)/h);%1*50




