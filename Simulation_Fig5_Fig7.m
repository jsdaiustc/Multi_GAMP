clear;
close all;
rand('seed',17)
randn('seed',17)
addpath(genpath(fileparts(mfilename('fullpath'))));
Params.random_seed = 10;             % The random state
Params.Fs            = 6400;         % The sampling frequency of the simulation signal
Params.N             = 64000;        % The length of the signal 
Params.Fn            =125;           % The fault characteristic frequency
Params.F             = 2000;         % The amplitude of the impulse
Params.mixture_ratio = [1,1];        % The mixing ratio of [impulses, harmonic, noise].
Params.noise_type = 'Laplacian';     % The noise type can be 'Gaussian' or 'Laplacian'
Fs=Params.Fs;

%% generate the 4-channel signals
window_size=[3200,4800,6400,9600];
TT=length(window_size);
C=zeros(9600,TT);
F2=zeros(9600,TT);
y_fft=zeros(9600,TT);
for tt=1:TT
window=window_size(tt);
Params.N= window;
[y,t,x] = Generate_Simulation_noseed(Params);
F2(1:window,tt)=[0:1:length(y)-1]'*Fs/length(y);
C(1:window,tt)=y;
f(:,tt)=y(1:3200);
y_envo= abs(hilbert(y));
y_fft(1:window,tt)= abs(fft(y_envo))/(length(y)/2);
y_h=  hilbert(y_envo);
yy(tt)=struct('cluster',y_h);
xx(tt)=struct('impsig',hilbert(abs(hilbert(x))));
end

%% plot the original signals 
figure(1);
[X,Y]=meshgrid(1:4,t);
Z=C;
plot3(X,Y,Z);
grid on;
xlabel('Channel');
ylabel('Time [s]');
zlabel('Amp.[m/s^2]');
title('a)')
set (gca,'position',[0.11,0.1,0.8,0.5],'FontSize',14,'FontName','Times New Rome' )



%% plot the envelope spectra of orignial signals
figure(2);
[X,Y]=meshgrid(1,F2(:,1));
Z=y_fft(:,1);
plot3(X,Y,Z);
grid on;hold on
[X,Y]=meshgrid(2,F2(:,2));
Z=y_fft(:,2);
plot3(X,Y,Z);
hold on;
[X,Y]=meshgrid(3,F2(:,3));
Z=y_fft(:,3);
plot3(X,Y,Z);
hold on;
[X,Y]=meshgrid(4,F2(:,4));
Z=y_fft(:,4);
plot3(X,Y,Z);
axis([0 4 50 800  0 0.3]);
xlabel('Channel');
ylabel('Frequency');
zlabel('Amp.[m/s^2]');
title('b)')
set (gca,'position',[0.12,0.1,0.8,0.5],'FontSize',14,'FontName','Times New Rome' )




figure(3);
fw=1000;             
%% GAMP
f_sample=[50:2:1000];  
[res_x,res_sample] =realMultichannel_GAMP(yy,f_sample,Fs,TT);
subplot(4,2,1);
line([125 250 375;125  250 375],[0 0 0;0.5 0.5 0.5],'color','k','linestyle','--'); hold on;
stem(res_sample,abs(res_x)/4,'marker','none','color','r');
axis([0 fw 0 0.5]);
ylabel('\fontname{Times New Roman}Amp.\fontname{Times New Roman}(m/s^2)');
title('a) RV-GAMP','fontname','Times New Roman');
set (gca,'position',[0.1,0.79,0.35,0.12],'FontSize',10,'FontName','Times New Rome')

%% MVMD
[u, u_hat, omega] = MVMD(C, 2000, 0, 3, 0, 1, 1e-7);
w=u(3,:,1)+u(2,:,1)+u(1,:,1);
N=length(C);
F = ([1:N]-1)*Fs/N;
F2= F(1:2001);
w_env=abs(fft(abs(hilbert(w))-mean(abs(hilbert(w)))))/(N/2);
subplot(4,2,2);
line([125 250 375;125  250 375],[0 0 0;0.5 0.5 0.5],'color','k','linestyle','--'); hold on;
plot(F2,  w_env(1:2001),'marker','none','color','blue')
axis([0 fw 0 0.2]);
title('b) MVMD','fontname','Times New Roman');
set (gca,'position',[0.52,0.79,0.35,0.12],'FontSize',10,'FontName','Times New Rome' )


%% MEMD 
imf = memd(C);
for nn=1:3200
   rr(nn,:)=imf(1,:,nn); 
end
rr=sum(rr,2);
emd_env=abs(fft(abs(hilbert(rr))-mean(abs(hilbert(rr)))))/(N/2);
subplot(4,2,3);
line([125 250 375;125  250 375],[0 0 0;0.5 0.5 0.5],'color','k','linestyle','--'); hold on;
plot(F2,  emd_env(1:2001),'marker','none','color','blue')
axis([0 fw 0 0.2]);
ylabel('\fontname{Times New Roman}Amp.\fontname{Times New Roman}(m/s^2)');
title('c) MEMD','fontname','Times New Roman');
set (gca,'position',[0.1,0.56,0.35,0.12],'FontSize',10,'FontName','Times New Rome' )


%% RV_ESPRIT
range=[0,1000];      
M=1000;              
[hat_f,~,hat_s]=RV_ESPRIT(y_h,M,Fs,range);
subplot(4,2,4);
line([125 250 375;125  250 375],[0 0 0;0.5 0.5 0.5],'color','k','linestyle','--'); hold on;
stem(hat_f,abs(hat_s),'marker','none','color','blue');
axis([0 fw 0 0.5]);
title('d) RV-ESPRIT','fontname','Times New Roman');
set (gca,'position',[0.52,0.56,0.35,0.12],'FontSize',10,'FontName','Times New Rome' )


%% SBL
[res_x1,res_sample1] =fault_frequency_learning(y_h,f_sample,Fs);
subplot(4,2,5);
line([125 250 375;125  250 375],[0 0 0;0.5 0.5 0.5],'color','k','linestyle','--'); hold on;
stem(res_sample1,abs(res_x1)/2,'marker','none','color','blue');
axis([0 fw 0 0.5]);
ylabel('\fontname{Times New Roman}Amp.\fontname{Times New Roman}(m/s^2)');
title('e) SBFL','fontname','Times New Roman');
set (gca,'position',[0.1,0.33,0.35,0.12],'FontSize',10,'FontName','Times New Rome' )


%% P-GSL
N=length(y);
F = ([1:N]-1)*Fs/N;
F2= F(1:2001);
[P_GSL_result] = P_GSL(y, Fs);
NMSE_GSL=(norm(x-P_GSL_result,'fro')/norm(x,'fro'))^2;
our_PSBL_enve=abs(fft(abs(hilbert(P_GSL_result)) -mean(abs(hilbert(P_GSL_result))) ))/(N/2);
subplot(4,2,6);
line([125 250 375;125  250 375],[0 0 0;0.5 0.5 0.5],'color','k','linestyle','--'); hold on;
plot(F2,  our_PSBL_enve(1:2001),'marker','none','color','blue')
axis([0 fw 0 0.5]);
title('f) P-GSL','fontname','Times New Roman');
set (gca,'position',[0.52,0.33,0.35,0.12],'FontSize',10,'FontName','Times New Rome' )


%% BPD
rho = 1;
K2 = round(N*0.04);
Method2.Name = 'WGL';
Method2.Initial_Size = 5;
Method2.SubName = 'MC';
Method2.gamma = 2;
Method2.window = 'gausswin';
Z2= IterGSS_improve(y, rho, K2, Method2);
y_BPD= real(Z2);
y_AdaESPGL_BPD=abs(fft(abs(hilbert(y_BPD)) -mean(abs(hilbert(y_BPD)))  ))/(N/2);
subplot(4,2,7);
line([125 250 375;125  250 375],[0 0 0;0.5 0.5 0.5],'color','k','linestyle','--'); hold on;
plot(F  ,y_AdaESPGL_BPD,'marker','none','color','blue')
axis([0 fw 0 0.5])
xlabel('\fontname{Times New Roman}Frequency\fontname{Times New Roman}(Hz)');
ylabel('\fontname{Times New Roman}Amp.\fontname{Times New Roman}(m/s^2)');
title('g) BPD')
set (gca,'position',[0.1,0.1,0.35,0.12],'FontSize',10,'FontName','Times New Rome' )



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
subplot(4,2,8);
line([125 250 375;125  250 375],[0 0 0;0.5 0.5 0.5],'color','k','linestyle','--'); hold on;
plot(F2,  y_AdaESPGL_enve(1:2001),'marker','none','color','blue')
axis([0 fw 0 0.5]);
xlabel('\fontname{Times New Roman}Frequency\fontname{Times New Roman}(Hz)');
title('h) AdaESPGL','fontname','Times New Roman');
set (gca,'position',[0.52,0.1,0.35,0.12],'FontSize',10,'FontName','Times New Rome' )







