function [cva,cvb,cvc,cvrho]=get_lscv(data,m,n,grid_length,hvec)
%test
load data100_n100_m100.mat
data=DATA{2};
hvec=0.0025:0.0025:0.06;
grid_length=50;m=100;n=m;

j=7
% percentage=100; 
%ngrid=100, 10% subset takes 1-2hour,est abc takes 4 min. 
%ngrid=100, 5% subset takes 45min.
% ngrid=50, 5% subset  takes 12min, 30s for one h, testing diff subsets
% ngrid=100, 5% subset almost get the same result as ngrid=50

% na_sample=100*percentage;
% randa_sample=sortrows(floor(rand(na_sample,2)*100+1));

[pair_cell1,pair_cell2]=meshgrid([1:n],[1:n]');
[pair_a1,pair_a2]=deal(reshape(pair_cell1,[],1),reshape(pair_cell2,[],1));
randa_sample=[pair_a1 pair_a2];
% nb_sample=10000*percentage;
nb_sample=1;
randb_sample=sortrows(floor(rand(nb_sample,3)*100+1));

T=1;
grid=[T/grid_length/2:T/grid_length:T];
tic
[paira,pairb,pairc,paira_same]=abc_pair_eval(data,m,n,grid_length,randa_sample,randb_sample);
toc

tic
for i=1:length(hvec)
    i
    h=hvec(i);
    int_u=cdf('Normal',(1-grid)/h,0,1)-cdf('Normal',(-grid)/h,0,1);
    edge=int_u'*int_u;
    %get esitmated t1-t5 and individual or daily t1-t5
%     tic
%     [A_leave,B_leave,C_leave,A,B,C]=est_abc_fun(data,m,n,h,grid_length,randa_sample,randb_sample);
        [A_leave,A]=est_a_fun(data,m,n,h,grid_length,randa_sample);
%     toc
    
%     tic
    %estimate cva, cvb, cvc
    cva1=sum(sum((A./edge/m/n/grid_length).^2));
%     cvb1=sum(sum((B./edge/m/n/(m-1)/grid_length).^2));
%     cvc1=sum(sum((C./edge/m/n/(n-1)/grid_length).^2));
    for k=1:na_sample
        estA_leave=A_leave(:,:,k)./edge/(m-1)/(n-1);
        cva2(k)=sum(estA_leave(paira{k}))-sum(estA_leave(paira_same{k}));
    end
%     for k=1:nb_sample
%         estB_leave=B_leave(:,:,k)./edge/(m-1)/(m-2)/(n-1);
%         cvb2(k)=sum(estB_leave(pairb{k}));
%     end
%     
%     for k=1:nb_sample
%         estC_leave=C_leave(:,:,k)./edge/(n-1)/(n-2)/(m-1);
%         cvc2(k)=sum(estC_leave(pairc{k}));
%     end
%     toc
    
    cva(i,j)=cva1-2*sum(cva2)/na_sample;
%     cvb(i,j)=cvb1-2*sum(cvb2)/nb_sample;
%     cvc(i,j)=cvc1-2*sum(cvc2)/nb_sample;   
end
toc


