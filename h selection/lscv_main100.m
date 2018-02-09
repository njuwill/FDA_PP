
% parpool('local100',50)
tic
n=100;m=n;
rep=300;
grid_length=100;
hvec=0.0025:0.0025:0.06;

tic
% load data100_n100_m100
load('C:\Users\mwang\Dropbox\finance project\pegasus\data\sim300\data300_nm100_rho5.mat')
for p=1:rep
    p
    data=DATA{p};
    [cva(p,:),cvb(p,:),cvc(p,:),cvrho(p,:)]=get_h(data,m,n,grid_length,hvec);
    
end
toc
save h_sim300_mn100 cva cvb cvc cvrho

clear all
tic
n=200;m=n;
rep=300;
grid_length=100;
hvec=0.0025:0.0025:0.06;

tic
% load data100_n100_m100
load('C:\Users\mwang\Dropbox\finance project\pegasus\data\sim300\data300_nm200_rho5.mat')
for p=1:rep
    data=DATA{p};
    [cva(p,:),cvb(p,:),cvc(p,:),cvrho(p,:)]=get_h(data,m,n,grid_length,hvec);
    
end
toc
save h_sim300_mn200 cva cvb cvc cvrho

clear all
tic
n=300;m=n;
rep=300;
grid_length=100;
hvec=0.0025:0.0025:0.06;

tic
% load data100_n100_m100
load('C:\Users\mwang\Dropbox\finance project\pegasus\data\sim300\data300_nm300_rho5.mat')
for p=1:rep
    p
    data=DATA{p};
    [cva(p,:),cvb(p,:),cvc(p,:),cvrho(p,:)]=get_h(data,m,n,grid_length,hvec);
    
end
toc
save h_sim300_mn300 cva cvb cvc cvrho