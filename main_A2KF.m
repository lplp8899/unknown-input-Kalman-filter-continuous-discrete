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

% this is the A2KF in the paper

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


%
[dim_sys,~]   =   size(x_real);
[dim_out,~]   =   size(y_real);
[dim_in,~] = size(u);
[dim_d,~]  =  size(d_real);
dim_sys_aug = dim_sys + dim_d;

%--------------  Augment the system matrices ------------
A_aug = [ A E;
        zeros(dim_d,dim_sys_aug-dim_d) zeros(dim_d)];
%
B_aug = [B; zeros(dim_d,dim_in);];
G_aug = [G zeros(dim_sys,dim_d); zeros(dim_d,dim_input_noise) eye(dim_d)]; % noise 
C_aug = [C zeros(dim_out,dim_d)];
%--------------  Augment the system matrices ------------

% initialization of vectors
inno=zeros(dim_out,gen);
x_ob_filter=zeros(dim_sys_aug,gen);
y_ob_filter=zeros(dim_out,gen);
P_ob_filter=zeros(dim_sys_aug,dim_sys_aug,gen);
residual=zeros(dim_out,gen);
K_gen = zeros(dim_sys_aug,dim_out,gen);


% for covariance estimation
S_inno=zeros(dim_out,dim_out,gen);
C_inno=zeros(dim_out,dim_out,gen);
norm_S_inno = zeros(1,gen); 
norm_C_inno = zeros(1,gen); 

Qk = Qk_real;
Rk = Rk_real;
Qdk = zeros(dim_d,dim_d);
Qdk_gen = zeros(dim_d,dim_d,gen);

% initial condition
x_ob_0 = [ zeros(dim_sys,1); zeros(dim_d,1)];
x_ob_0 = [ 10; 10; 10; 10; zeros(dim_d,1)];
P_x0 = 1e2*eye(dim_sys_aug);


% start of the filter
for k=1:gen/1

    % predict
%     x_ob=x_ob_0+model_fx(x_ob_0,u(:,k),g(:,k))*delta_t;
	x_ob = x_ob_0 + (A_aug*x_ob_0 + B_aug*u(:,k) )*dt;
   
    % Adaptive covariance
%         Qk_aug=[ Qk_real, zeros(dim_sys,dim_d); zeros(dim_d,dim_sys), diag([Qdk(1,1), Qdk(2,2) ])];
        Qk_aug=[ Qk_real, zeros(dim_sys,dim_d); zeros(dim_d,dim_sys), diag( diag(Qdk) )];
    
%     P = A_aug*P_x0*A_aug' + G_aug*Qk_aug*G_aug';
    Ad_aug =  eye(dim_sys_aug) + A_aug*dt;
    % prediction error covariance
    P = Ad_aug*P_x0*Ad_aug' + G_aug*Qk_aug*G_aug'*dt;
       
    V=C_aug*P*C_aug'+Rk;

    % Kalman gain
    K=P*C_aug'/V;
    K_gen(:,:,k) = K;

%     z_ob=model_output(x_ob);
    z_ob = C_aug*x_ob;
    % innovation
    inno(:,k)=y_real(:,k)-z_ob;

    % estimation    
    x_ob_filter(:,k)=x_ob+K*inno(:,k);
    P_ob_filter(:,:,k)=(eye(dim_sys_aug)-K*C_aug)*P*(eye(dim_sys_aug)-K*C_aug)'+K*Rk*K';
%     tr_P(k)=trace(P_ob_filter(:,:,k));
    
    %next real input
    x_ob_0=x_ob_filter(:,k);
    P_x0=P_ob_filter(:,:,k);


    %------------------ adaptive Qd estimation start ------------------%
    S_inno(:,:,k) = inno(:,k)*inno(:,k)';
    
    N_win = 1;% width of moving window
    % compute innovation covariance
    if N_win == 1
        C_inno(:,:,k) = S_inno(:,:,k);
    end
    if N_win > 1
        if k == 1
            C_inno(:,:,k) = S_inno(:,:,k);
        elseif k >1 && k < N_win
            C_inno(:,:,k) = sum(S_inno(:,:,1:k),3)/(k-1);
        elseif k > N_win
            C_inno(:,:,k) = sum(S_inno(:,:,k-N_win+1:k),3)/(N_win-1);
        end 
    end

    % compute the Qd
    Qf_0 = ( C_inno(:,:,k) - C*G*Qk*G'*C'*dt - Rk );
    Qf_part = diag( abs(diag(Qf_0)) ); % only use the diagonal elements of Qf_0
    CE_pinv = pinv(C*E);

    % estimation of Qd
    Qdk = CE_pinv* Qf_part * CE_pinv'/dt/dt; % Eq(60) in the paper, Ed=E*dt

    %--------- the following can also compute Qdk
    % or use Ed (the discretized version of E)
%     Ed = E*dt;
%     CEd_pinv = pinv(C*Ed);
%     Qdk = CEd_pinv* Qf_part * CEd_pinv';
    %------------------ adaptive Qd estimation end ------------------%
    Qdk_gen(:,:,k) = Qdk;


end



% results
% 
Time=dt*(1:k);

% Qd
figure;
subplot(211);hold on;plot(Time,squeeze(Qdk_gen(1,1,1:k)),'b','linewidth',2); ylabel('Qd_{11} ','fontsize',12);grid;set(gca,'fontsize',12); %plot(Time,P_froot(1,1:k));plot(Time,-P_froot(1,1:k));
set(gca,'ylim',[0 1]);% title('Error of Fault Estimation')
if case_no == 3
    set(gca,'ylim',[0 5]);
end
subplot(212);hold on;plot(Time,squeeze(Qdk_gen(2,2,1:k)),'b','linewidth',2); ylabel('Qd_{22}','fontsize',12);grid;set(gca,'fontsize',12); %plot(Time,P_froot(2,1:k));plot(Time,-P_froot(2,1:k));
set(gca,'ylim',[0 0.1]);
if case_no == 2
    set(gca,'ylim',[0 0.4]);
end
if case_no == 3
    set(gca,'ylim',[0 1.5]);
end
h1=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time (s)','fontsize',12)
set(h1,'Box','off')
% end 

% x
figure;
subplot(221); hold on; plot(Time,x_real(1,1:k),'r','linewidth',2); plot(Time,x_ob_filter(1,1:k),'b--','linewidth',2); ylabel('x_1','fontsize',12);grid;set(gca,'fontsize',12);
subplot(222); hold on; plot(Time,x_real(2,1:k),'r','linewidth',2); plot(Time,x_ob_filter(2,1:k),'b--','linewidth',2); ylabel('x_2','fontsize',12);grid;set(gca,'fontsize',12);
subplot(223); hold on; plot(Time,x_real(3,1:k),'r','linewidth',2); plot(Time,x_ob_filter(3,1:k),'b--','linewidth',2); ylabel('x_3','fontsize',12);grid;set(gca,'fontsize',12);
subplot(224); hold on; plot(Time,x_real(4,1:k),'r','linewidth',2); plot(Time,x_ob_filter(4,1:k),'b--','linewidth',2); ylabel('x_4','fontsize',12);grid;set(gca,'fontsize',12);
h=legend('True','Estimation'); set(h,'color','none','edgecolor','white');
% h1=axes('position',[0.10 0.05 0.84 0.86],'fontsize',12);axis off;title(h1,'Fault Free')
h2=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time (s)','fontsize',12)
set(h2,'Box','off')

% P
for i=1:dim_sys
    P_sroot(i,:)=sqrt(P_ob_filter(i,i,:));
end
% error of estimation
figure;
subplot(221);hold on;plot(Time,x_real(1,1:k)-x_ob_filter(1,1:k),'b','linewidth',2); ylabel('\Deltax_1 ','fontsize',12); grid;set(gca,'fontsize',12);%plot(Time,P_sroot(1,1:k));plot(Time,-P_sroot(1,1:k)); axis([1 k/100 -9e-2 9e-2])
subplot(222);hold on;plot(Time,x_real(2,1:k)-x_ob_filter(2,1:k),'b','linewidth',2); ylabel('\Deltax_2 ','fontsize',12); grid;set(gca,'fontsize',12);%plot(Time,P_sroot(2,1:k));plot(Time,-P_sroot(2,1:k)); axis([1 k/100 -1e-3 1e-3])
subplot(223);hold on;plot(Time,x_real(3,1:k)-x_ob_filter(3,1:k),'b','linewidth',2); ylabel('\Deltax_3 ','fontsize',12); grid;set(gca,'fontsize',12);%plot(Time,P_sroot(2,1:k));plot(Time,-P_sroot(2,1:k)); axis([1 k/100 -1e-3 1e-3])
subplot(224);hold on;plot(Time,x_real(4,1:k)-x_ob_filter(4,1:k),'b','linewidth',2); ylabel('\Deltax_4 ','fontsize',12); grid;set(gca,'fontsize',12);%plot(Time,P_sroot(2,1:k));plot(Time,-P_sroot(2,1:k)); axis([1 k/100 -1e-3 1e-3])
h1=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time (s)','fontsize',12)
set(h1,'Box','off')

% RMSE
% compute the RMSE after 10 steps. 
% so that the initialization error will not influence the RMSE
    fprintf(1,'the RMSE of state estimation using AEKF is:\n%f\n%f\n%f\n%f\n ',...
        sqrt(mean((x_real(1,10:k)-x_ob_filter(1,10:k)).^2)),...
        sqrt(mean((x_real(2,10:k)-x_ob_filter(2,10:k)).^2)),...
        sqrt(mean((x_real(3,10:k)-x_ob_filter(3,10:k)).^2)),...
        sqrt(mean((x_real(4,10:k)-x_ob_filter(4,10:k)).^2)));
    % fprintf(1,'the norm of root mean square error is:\n%f\n ',norm(RMSE(:)))

    fprintf(1,'the RMSE of unknown input estimation using AEKF is:\n%f\n%f\n%f\n ',...
        sqrt(mean((d_real(1,10:k)-x_ob_filter(dim_sys+1,10:k)).^2)),...
        sqrt(mean((d_real(2,10:k)-x_ob_filter(dim_sys+2,10:k)).^2)));

% if dim_sys > 6
% x-augmented states, i.e. the unknown inputs
figure;
subplot(211);hold on; plot(Time,d_real(1,1:k),'r','linewidth',2);plot(Time,x_ob_filter(dim_sys+1,1:k),'b--','linewidth',2);  ylabel('d_1 ','fontsize',12);grid;set(gca,'fontsize',12);
set(gca,'ylim',[-0.2 0.8]);% title('Fault Estimation')
subplot(212);hold on; plot(Time,d_real(2,1:k),'r','linewidth',2);plot(Time,x_ob_filter(dim_sys+2,1:k),'b--','linewidth',2);  ylabel('d_2','fontsize',12);grid;set(gca,'fontsize',12);
set(gca,'ylim',[-0.45 0.45]);
h2=legend('True','Estimation'); set(h2,'color','white','edgecolor','black');
% h1=axes('position',[0.10 0.05 0.84 0.86],'fontsize',12);axis off;title(h1,'Fault ')
h1=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time (s)','fontsize',12)
set(h1,'Box','off')
    

