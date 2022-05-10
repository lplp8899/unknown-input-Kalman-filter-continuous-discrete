%----------------------------------------------------------------------
% This is the code for the following paper:
%
%      P. Lu. Comparative Study of State and Unknown Input Estimation for
%      Continuous-discrete Stochastic Systems. Franklin Open, 2022, 1-14.
%
%   This program shows the results of Section 6 of the paper.
%
%   Author: Peng Lu, the University of Hong Kong
%
%   email:  lupeng@hku.hk
%
%   released in April 2022
%
%   Please cite the above paper if you find this code useful
%----------------------------------------------------------------------

% this is the R4SKF in the paper

% x= Ax+Bu+Ef
% y= Hx

close all
clear all
% clc

% choose different simulation cases:
% 1 for normal case, 2 for case 2 and 3 for case 3,
% if you failed to input the correct number, the default case is case 1.
disp('1. case 1   2.case 2: faster dynamics   3.case 3: large measurement noise');
case_no = input('Choose 1, 2 or 3: ');

% generate system data including measurements
% real states are used to compute the estimation errors
gene_sys_data
      

%------------------------------ initialization  ------------------------------  
% x_est_0 = zeros(dim_sys,1); % the true initial condition, which is not used
x_est_0 = [10; 10; 10; 10];
P_est_0=1e0*eye(dim_sys);

Qk = Qk_real;
Rk = Rk_real;

%-------------------------------------------------------------------------

%
x_pred = zeros(dim_sys,1);
inno = zeros(dim_out,gen); % innovation vector

error_x = zeros(dim_sys,gen);

x_est = zeros(dim_sys,gen);
y_est = zeros(dim_out,gen);
P_est = zeros(dim_sys,dim_sys,gen);
residual = zeros(dim_out,gen);

d_est=zeros(dim_d,gen);
Pd_est=zeros(dim_d,dim_d,gen);



for k=1:gen/1

    %-------------------- step 1: prediction without UI
   
    % obtain discrete matrices
    Ad = eye(dim_sys)+A*dt;
    Bd = B*dt;
    Ed = E*dt;
    Gd = G*dt;
    x_pred_noUI = Ad*x_est_0 + Bd*u(:,k);
    % you can also use the following to compute x_pred_noUI:
    x_pred_noUI = x_est_0 + ( A*x_est_0 + B*u(:,k) )*dt;
    
    P_pred = Ad*P_est_0*Ad' + G*Qk*G'*dt;

    
    %-------------------- step 2: UI estimation
    inno_noUI = y_real(:,k) - C*x_pred_noUI;
    CEk = C*Ed;
    R_hat = C*P_pred*C'+Rk;
    inv_R = eye(dim_out)/R_hat;
    Fd = inv(CEk'*CEk)*CEk'; % get the MP inverse
    d_est(:,k) = Fd*inno_noUI;
    Pd_est = inv(CEk'*inv_R*CEk);
    
    
    %-------------------- step 3: state estimation
    x_pred_UI = x_est_0 + ( A*x_est_0 + B*u(:,k) + E*d_est(:,k) )*dt;


    %-------------------- step 4: measurement update
    % use standard Kalman gain
    Kk = P_pred*C'*inv_R;
 
    Lk = Kk+(eye(dim_sys)-Kk*C)*Ed*Fd;
    % innovation
    inno(:,k) = y_real(:,k) - C*x_pred_UI;
    % final state estimation
    x_est(:,k) = x_pred_UI + Kk * inno(:,k);
%     % you can also use this one
%     x_est(:,k) = x_pred_noUI + Lk*inno_noUI;
    % covariance
    P_est(:,:,k) = (eye(dim_sys)-Lk*C)*P_pred*(eye(dim_sys)-Lk*C)'+Lk*Rk*Lk';   
    
    % next iteration
    x_est_0 = x_est(:,k);
    P_est_0 = P_est(:,:,k);

    
end


% 
Time=dt*(1:k);

% x
figure;
subplot(411); hold on; plot(x_real(1,1:k),'r','linewidth',2); plot(x_est(1,1:k),'b--','linewidth',2); ylabel('x_1','fontsize',12);grid; set(gca,'fontsize',12);
subplot(412); hold on; plot(x_real(2,1:k),'r','linewidth',2); plot(x_est(2,1:k),'b--','linewidth',2); ylabel('x_2','fontsize',12);grid;set(gca,'fontsize',12); 
subplot(413); hold on; plot(x_real(3,1:k),'r','linewidth',2); plot(x_est(3,1:k),'b--','linewidth',2); ylabel('x_3','fontsize',12);grid;set(gca,'fontsize',12);
subplot(414); hold on; plot(x_real(4,1:k),'r','linewidth',2); plot(x_est(4,1:k),'b--','linewidth',2); ylabel('x_4','fontsize',12);grid;set(gca,'fontsize',12);
h=legend('True','Estimation'); set(h,'color','white','edgecolor','black');
% h1=axes('position',[0.10 0.05 0.84 0.86],'fontsize',12);axis off;title(h1,'Fault Free')
h2=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time','fontsize',12)
set(h2,'Box','off')

% P
% error of estimation
figure;
subplot(411);hold on;plot(Time,x_real(1,1:k)-x_est(1,1:k),'b','linewidth',2); ylabel('\Deltax_1 ','fontsize',12); grid;set(gca,'fontsize',12);%plot(Time,P_sroot(1,1:k));plot(Time,-P_sroot(1,1:k)); axis([1 k/100 -9e-2 9e-2])
% set(gca,'ylim',[-0.01 0.01]);
subplot(412);hold on;plot(Time,x_real(2,1:k)-x_est(2,1:k),'b','linewidth',2); ylabel('\Deltax_2 ','fontsize',12); grid;set(gca,'fontsize',12);%plot(Time,P_sroot(2,1:k));plot(Time,-P_sroot(2,1:k)); axis([1 k/100 -1e-3 1e-3])
% set(gca,'ylim',[-0.01 0.01]);    
subplot(413);hold on;plot(Time,x_real(3,1:k)-x_est(3,1:k),'b','linewidth',2); ylabel('\Deltax_3 ','fontsize',12); grid;set(gca,'fontsize',12);%plot(Time,P_sroot(2,1:k));plot(Time,-P_sroot(2,1:k)); axis([1 k/100 -1e-3 1e-3])
subplot(414);hold on;plot(Time,x_real(4,1:k)-x_est(4,1:k),'b','linewidth',2); ylabel('\Deltax_4 ','fontsize',12); grid;set(gca,'fontsize',12);%plot(Time,P_sroot(2,1:k));plot(Time,-P_sroot(2,1:k)); axis([1 k/100 -1e-3 1e-3])
h1=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time','fontsize',12)
set(h1,'Box','off')

% RMSE root mean square error 
% compute the RMSE after 10 steps. 
% so that the initialization error will not influence the RMSE
fprintf(1,'the RMSE of state estimation using R4KF is:\n%f\n%f\n%f\n%f\n ',...
    sqrt(mean((x_real(1,10:k)-x_est(1,10:k)).^2)),...
    sqrt(mean((x_real(2,10:k)-x_est(2,10:k)).^2)),...
    sqrt(mean((x_real(3,10:k)-x_est(3,10:k)).^2)),...
    sqrt(mean((x_real(4,10:k)-x_est(4,10:k)).^2)));

fprintf(1,'the RMSE of unknown input estimation using R4SKF is:\n%f\n%f\n%f\n ',...
    sqrt(mean((d_real(1,10:k)-d_est(1,10:k)).^2)),...
    sqrt(mean((d_real(2,10:k)-d_est(2,10:k)).^2)));

% unknown input estimation
figure;
subplot(211);hold on; plot(Time,d_real(1,1:k),'r','linewidth',2);plot(Time,d_est(1,1:k),'b--','linewidth',2);  ylabel('d_1 ','fontsize',12);grid;set(gca,'fontsize',12);
%     set(gca,'ylim',[-0.3 0.8]);% title('Fault Estimation')
set(gca,'ylim',[-3 3]);% title('Fault Estimation')
subplot(212);hold on; plot(Time,d_real(2,1:k),'r','linewidth',2);plot(Time,d_est(2,1:k),'b--','linewidth',2);  ylabel('d_2','fontsize',12);grid;set(gca,'fontsize',12);
%     set(gca,'ylim',[-0.5 0.5]);
set(gca,'ylim',[-2 2]);
h2=legend('True','Estimation'); set(h2,'color','white','edgecolor','black');
% h1=axes('position',[0.10 0.05 0.84 0.86],'fontsize',12);axis off;title(h1,'Fault ')
h1=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time','fontsize',12)
set(h1,'Box','off')
    



