clear;
close all;
addpath(genpath(fileparts(mfilename('fullpath'))));
%% IMS bearing dataset
y1=load('2004.02.16.06.02.39_550');
y1=y1(:,1);
y2=load('2004.04.16.19.42.55_6140');
y2=y2(:,3);
y_sum=[y1,y2];
Fs=20480;

%% create the 2-channel vibration signals
aa1=3750; aa2=5120;
window_size=[aa1,aa2];
TT=length(window_size);
C=zeros(aa2,TT);
F2=zeros(aa2,TT);
y_fft=zeros(aa2,TT);
for tt=1:TT
    window=window_size(tt);
    Params.N= window;
    y=y_sum(:,tt);
    y=y(1:window);
    F2(1:window,tt)=[0:1:length(y)-1]'*Fs/length(y);
    C(1:window,tt)=y;
    y_envo= abs(hilbert(y));
    y_fft(1:window,tt)= abs(fft(y_envo))/(length(y)/2);
    y_h=  hilbert(y_envo);
    yy(tt)=struct('cluster',y_h);
end
N=length(y);
t = (0 : N-1) / Fs;


%% generate the 2-channel signals
figure(1);
[X,Y]=meshgrid(1:2,t);
Z=C;
plot3(X,Y,Z);
grid on;
axis([0 2 0 0.25]);
xlabel('Channel');
ylabel('Time [s]');
zlabel('Amp.[m/s^2]');
title('a)')
set (gca,'position',[0.14,0.1,0.8,0.5],'FontSize',14,'FontName','Times New Rome' )


%% plot the envelope spectra of orignial signals
figure(2);
[X,Y]=meshgrid(1,F2(:,1));
Z=y_fft(:,1);
plot3(X,Y,Z);
grid on;hold on
[X,Y]=meshgrid(2,F2(:,2));
Z=y_fft(:,2);
plot3(X,Y,Z);
axis([0 2 50 1850  0 0.02]);
xlabel('Channel');
ylabel('Frequency [Hz]');
zlabel('Amp.[m/s^2]');
title('b)')
set (gca,'position',[0.15,0.2,0.8,0.5],'FontSize',14,'FontName','Times New Rome' )

figure(3);
fw=1850;             %%%%%%%坐标横轴大小
%% GAMP
f_sample=[50:2:1000];    
[res_x,res_sample] =realMultichannel_GAMP(yy,f_sample,Fs,TT);
subplot(3,2,1);
line([236.4 472.8 709.2;236.4 472.8 709.2],[0 0 0;1.5 1.5 1.5],'color','k','linestyle','--'); hold on;
stem(res_sample,abs(res_x)/4,'marker','none','color','r');
axis([0 fw 0 0.02]);
ylabel('\fontname{Times New Roman}Amp.\fontname{Times New Roman}(m/s^2)');
title('a) RV-GAMP','fontname','Times New Roman');
set (gca,'position',[0.1,0.65,0.35,0.15],'FontSize',10,'FontName','Times New Rome')

%% MVMD
[u, u_hat, omega] = MVMD(C, 2000, 0, 3, 0, 1, 1e-7);
w=u(3,:,1)+u(2,:,1)+u(1,:,1);
N=length(C);
F = ([1:N]-1)*Fs/N;
F2= F(1:2001);
w_env=abs(fft(abs(hilbert(w))-mean(abs(hilbert(w)))))/(N/2);
subplot(3,2,2);
line([236.4 472.8 709.2;236.4 472.8 709.2],[0 0 0;1.5 1.5 1.5],'color','k','linestyle','--'); hold on;
plot(F2,  w_env(1:2001),'marker','none','color','blue')
axis([0 fw 0 0.02]);
title('b) MVMD','fontname','Times New Roman');
set (gca,'position',[0.52,0.65,0.35,0.15],'FontSize',10,'FontName','Times New Rome' )


%% RV_ESPRIT
range=[0,1000];      %%%%%%滤波范围 
M=1000;               %%%%%%窗口大小
[hat_f,~,hat_s]=RV_ESPRIT(y_h,M,Fs,range);
subplot(3,2,3);
line([236.4 472.8 709.2;236.4 472.8 709.2],[0 0 0;1.5 1.5 1.5],'color','k','linestyle','--'); hold on;
stem(hat_f,abs(hat_s),'marker','none','color','blue');
axis([0 fw 0 0.02]);
ylabel('\fontname{Times New Roman}Amp.\fontname{Times New Roman}(m/s^2)');
title('c) RV-ESPRIT','fontname','Times New Roman');
set (gca,'position',[0.1,0.4,0.35,0.15],'FontSize',10,'FontName','Times New Rome' )


%% SBL
[res_x1,res_sample1] =fault_frequency_learning(y_h,f_sample,Fs);
subplot(3,2,4);
line([236.4 472.8 709.2;236.4 472.8 709.2],[0 0 0;1.5 1.5 1.5],'color','k','linestyle','--'); hold on;
stem(res_sample1,abs(res_x1)/2,'marker','none','color','blue');
axis([0 fw 0 0.02]);
title('d) SBFL','fontname','Times New Roman');
set (gca,'position',[0.52,0.4,0.35,0.15],'FontSize',10,'FontName','Times New Rome' )


%% P-GSL
N=length(y);
F = ([1:N]-1)*Fs/N;
F2= F(1:2001);
[P_GSL_result] = P_GSL(y, Fs);
our_PSBL_enve=abs(fft(abs(hilbert(P_GSL_result)) -mean(abs(hilbert(P_GSL_result))) ))/(N/2);
subplot(3,2,5);
line([236.4 472.8 709.2;236.4 472.8 709.2],[0 0 0;1.5 1.5 1.5],'color','k','linestyle','--'); hold on;
plot(F2,  our_PSBL_enve(1:2001),'marker','none','color','blue')
axis([0 fw 0 0.02]);
 xlabel('\fontname{Times New Roman}Frequency\fontname{Times New Roman}(Hz)');
 ylabel('\fontname{Times New Roman}Amp.\fontname{Times New Roman}(m/s^2)');
title('e) P-GSL','fontname','Times New Roman');
set (gca,'position',[0.1,0.15,0.35,0.15],'FontSize',10,'FontName','Times New Rome' )


%% AdaESPGL, downloaded from https://zhaozhibin.github.io/  %%%%%%%%%%%%%%%%%%%%%%%%%
Params.Fs            = Fs;     % The sampling frequency of the simulation signal
Params.N             = N;      % The length of the signal
Params.N1    = 4;              % The samples of one impulse
Params.M     = 4;              % The number of periods
Params.Fn_N  = 0;              % a vector which contains the period of each component (Fs / fc)
Params.mu    = 9.235e-4;       % The parameter related to sparsity within groups
Params.pen   = 'atan';         % The penalty function
Params.rho   = 1;              % The degree of nonconvex
Params.Nit   = 100;            % The number of iteration
% Estimate noise
[C,L]=wavedec(y,5,'sym8');
c1=detcoef(C,L,1);
est_noise=median(abs(c1-median(c1)))/0.678;
Params.lam= 0.272*est_noise + 0.044;
[AdaESPGL_result] = AdaESPGL(y, Params);
y_AdaESPGL_enve= abs(fft(abs(hilbert(AdaESPGL_result))-mean(abs(hilbert(AdaESPGL_result)))))/(N/2);
subplot(3,2,6);
line([59.6 119.2 178.8;59.6 119.2 178.8],[0 0 0;0 0 0],'color','k','linestyle','--'); hold on;
plot(F2,  y_AdaESPGL_enve(1:2001),'marker','none','color','blue')
text(800, 0.01, 'Failed');
axis([0 fw 0 0.02]);
xlabel('\fontname{Times New Roman}Frequency\fontname{Times New Roman}(Hz)');
title('f) AdaESPGL','fontname','Times New Roman');
set (gca,'position',[0.52,0.15,0.35,0.15],'FontSize',10,'FontName','Times New Rome' )







