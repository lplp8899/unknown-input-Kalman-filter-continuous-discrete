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

% continuous-discrete stochastic sys model:
% x_dot(t) = A(t)*x(t) + B(t)*u(t) + E(t)*d(t) + G(t)*w(t);
% y(k+1) = C(k+1)*x(k+1)  + v(k+1)

% close all
% clear all
% clc


gen = 1000; % number of iterations


% ------------------ one example
%  unstable
A=[ 1.9527  -0.0075 0.0663  0.0437;
    0.0017  1.0452  0.0056  -0.0242;
    0.0092  0.0064 -0.1975 0.00128;
    0 0 1 0;];

B=[ 0.554 0.156;
    0.246 -0.982;
    0.320 0.560;
    0 0;];

E = B; % disturbance distribution matrix
% E = eye(4);

G = eye(size(A,1));

C=[1 0 0 0;
    0 1 0 0;
    0 0 0 1];
% C = eye(4);

D = zeros(size(C,1),size(B,2));

% disp('If augment or TEKF, eig(A) should be smaller than 1 for discrete system');
% disp('If augment or TEKF, eig(A) should have negative real parts for continuous system (Herwitz)');
% 
% eig_A = eig(A)
% if find(eig_A >= 0) % did not check complex numbers
%     disp('NOT a stable system matrix');
% end 

% dimensions
dim_sys = size(A,1);
dim_out = size(C,1);
dim_d = size(E,2);
dim_in = size(B,2);
dim_input_noise = size(G,2);

x_real=zeros(dim_sys,gen);
y_real=zeros(dim_out,gen);

% process and measurement noise covariance
Qk_real = diag([1e-6, 1e-6, 1e-6, 1e-6]);

% normal measurement noise
% Rk_real = diag([1e-6, 1e-6, 1e-6, 1e-6  ]) *1e-1;
Rk_real = diag([1e-6, 1e-6, 1e-6  ]) *1e-1;
% case 3: Larger measurement noise
if case_no == 3
    Rk_real = diag([1e-6, 1e-6, 1e-6  ]) *1e1;
end


dt = 0.01;
% input is the random input
randn('state',1000)
w = sqrt(Qk_real)*randn(dim_input_noise,gen);
randn('state',2000)
v = sqrt(Rk_real)*randn(dim_out,gen);


% real unknown inputs
dy = 0.001;
d_real = zeros(dim_d,gen);
d_real(1,1:200) = 0;
d_real(1,301:700) = 0.5;
freq = 0.5;
% case 2: faster unknown input dynamics
if case_no == 2
    freq = 5;
end
d_real(2,201:600)=d_real(2,201:600)+ 0.4*sin(2*pi*freq*0.01*(1:400));

% input
u = zeros(dim_in,gen);
% u = ones(dim_in,gen);
% u(1,:) = 0.0;
% u(2,:) = 0.000;
% u(3,:) = -9.8065;
dim_in = size(u,1);

% real initial state
x_real_0 = zeros(dim_sys,1);
x_init_real = x_real_0;

% generate real state and measurement data
for k=1:gen
    % continuous systems
    x_real_dot = A*x_real_0 + B*u(:,k) + E*d_real(:,k) + G*w(:,k);
    x_real(:,k) = x_real_0 + x_real_dot*dt;
    y_real(:,k) = C*x_real(:,k) + v(:,k);
    x_real_0 = x_real(:,k);
end
Time = dt*(1:k);



%--------- observability check
Ob = obsv(A,C);
unobs = rank(A) - rank(Ob) % number of unobservable states




% % plots
% figure;
% subplot(411); hold on; grid; plot(Time,x_real(1,1:k),'b','linewidth',1.5); 
% % ylabel('f_V (m/s)','fontsize',13);grid; set(gca,'ylim',[-2.5 2.5],'fontsize',13)
% subplot(412); hold on; grid; plot(Time,x_real(2,1:k),'b','linewidth',1.5); 
% % ylabel('f_{\beta} (rad)','fontsize',13);grid;set(gca,'ylim',[-0.1 0.1],'fontsize',13)
% subplot(413); hold on; grid; plot(Time,x_real(3,1:k),'b','linewidth',1.5); 
% % ylabel('f_{\beta} (rad)','fontsize',13);grid;set(gca,'ylim',[-0.1 0.1],'fontsize',13)
% subplot(414); hold on; grid; plot(Time,x_real(4,1:k),'b','linewidth',1.5); 
% % ylabel('f_{\beta} (rad)','fontsize',13);grid;set(gca,'ylim',[-0.1 0.1],'fontsize',13)
% % legend('estimation','true','fontsize',13);
% % h1=axes('position',[0.50 0 0.0001 0.0001],'fontsize',13);axis off;
% % title('     time (s)','fontsize',13)
% 
% 
% figure;
% subplot(411); hold on; grid; plot(Time,y_real(1,1:k),'b','linewidth',1.5); plot(Time,x_real(1,1:k),'r','linewidth',1.5);
% % ylabel('f_V (m/s)','fontsize',13);grid; set(gca,'ylim',[-2.5 2.5],'fontsize',13)
% subplot(412); hold on; grid; plot(Time,y_real(2,1:k),'b','linewidth',1.5); plot(Time,x_real(2,1:k),'r','linewidth',1.5);
% % ylabel('f_{\beta} (rad)','fontsize',13);grid;set(gca,'ylim',[-0.1 0.1],'fontsize',13)
% subplot(413); hold on; grid; plot(Time,y_real(3,1:k),'b','linewidth',1.5); plot(Time,x_real(3,1:k),'r','linewidth',1.5);
% % ylabel('f_{\beta} (rad)','fontsize',13);grid;set(gca,'ylim',[-0.1 0.1],'fontsize',13)
% subplot(414); hold on; grid; plot(Time,y_real(4,1:k),'b','linewidth',1.5); plot(Time,x_real(4,1:k),'r','linewidth',1.5);
% % ylabel('f_{\beta} (rad)','fontsize',13);grid;set(gca,'ylim',[-0.1 0.1],'fontsize',13)
% legend('measurement','state','fontsize',13);
% % h1=axes('position',[0.50 0 0.0001 0.0001],'fontsize',13);axis off;
% % title('     time (s)','fontsize',13)

% dk
% figure;
% subplot(211);hold on;plot(Time,d_real(1,1:k),'b','linewidth',2); ylabel('d_1 ','fontsize',12);grid;set(gca,'fontsize',12); %plot(Time,P_froot(1,1:k));plot(Time,-P_froot(1,1:k));
% set(gca,'ylim',[-0.2 0.6]);
% subplot(212);hold on;plot(Time,d_real(2,1:k),'b','linewidth',2); ylabel('d_2','fontsize',12);grid;set(gca,'fontsize',12); %plot(Time,P_froot(2,1:k));plot(Time,-P_froot(2,1:k));
% % set(gca,'ylim',[-0.5 0.5]);
% h1=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time','fontsize',12)
% set(h1,'Box','off')





