clear all
load xiy_logis_sim100_mn100
load xix_logis_sim300_mn200 
load xiy_logis_sim300_mn200 
load('C:\Users\mwang\Dropbox\finance project\pegasus\data\sim300\data300_nm100_rho5.mat')

%% compare all the esitmated and true xi with scatter plot
rep=300;
xi_y_est=[];
xi_x_est=[];
for i=1:rep
   xi_x_tmp=xi_xb_est{i}; 
   xi_y_tmp=xi_yb_est{i};
   xi_y_est=[xi_y_est; xi_y_tmp];
   xi_x_est=[xi_x_est; xi_x_tmp];
end

xi_x_true=[];
xi_y_true=[];
for i=1:rep
   xi_true_tmp=XI{i};
   xi_x_tmp=[xi_true_tmp{1}];
   xi_y_tmp=[xi_true_tmp{2}' xi_true_tmp{3}'];
   xi_x_true=[xi_x_true; xi_x_tmp];
   xi_y_true=[xi_y_true; xi_y_tmp];
end

subplot(1,3,1);scatter(xi_x_est(:,1), xi_x_true(:,1),'.');xlim([-4,4]);ylim([-4,4]);xlabel('$\widehat{\xi_{i1}^X}$','Interpreter','latex');ylabel('$\xi_{i1}^X$','Interpreter','latex')
subplot(1,3,2);scatter(xi_y_est(:,1), xi_y_true(:,1),'.');xlim([-4,4]);ylim([-4,4]);xlabel('$\widehat{\xi_{j1}^Y}$','Interpreter','latex');ylabel('$\xi_{j1}^Y$','Interpreter','latex')
subplot(1,3,3);scatter(xi_y_est(:,2), xi_y_true(:,2),'.');xlim([-4,4]);ylim([-4,4]);xlabel('$\widehat{\xi_{j2}^Y}$','Interpreter','latex');ylabel('$\xi_{j2}^Y$','Interpreter','latex')

%estimate the mean and var first, then compare to the true value
for i=1:rep
    x_mean(i,:)=mean(parax{i});
    x_var(i,:)=var(parax{i});
end

%test
xi_test=paray_logistic{1};
xi_test1=paray_logistic1{1};

%% get the diffenence of true and esitimated xi


difflag_xiy=[];
for i=1:rep
    xi_y=XI{i};
    xi_true=[xi_y{2}; xi_y{3}]';
    xi_ylag=paray_logistic1{i};
    xi_tmp=xi_ylag(:,1:2)-xi_true;
    difflag_xiy=[difflag_xiy;xi_tmp];
%     ylag_mean(i,:)=mean(xi_ylag(:,1:2)-xi_true);
%     ylag_var(i,:)=var(xi_ylag(:,1:2)-xi_true);
end

[mean(difflag_xiy);mean(diff_xiy)]
[var(difflag_xiy);var(diff_xiy)]


for i=1:m
    ylag_mean(i,:)=mean(paray_logistic1{i});
    ylag_var(i,:)=var(paray_logistic1{i});
end
% [mean(ylag_mean);mean(ylag_var)]
[mean(y_mean);mean(y_var)]

% [var(ylag_mean);var(ylag_var)]
% [var(y_mean);var(y_var)]


% xi_test=parax_logistic{1};
xi_sim=XI{1};
xi_true=xi_sim{3};
scatter(xi_test1(:,2),xi_true)
[xi_test1(:,2) xi_true']

scatter(xi_test(:,1),xi_test1(:,1))
scatter(xi_test(:,1),xi_origianl(:,1))