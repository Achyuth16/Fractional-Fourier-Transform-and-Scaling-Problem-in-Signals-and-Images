clc;clear all;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=100;
t=(0:1/fs:1-(1/fs));
a=1.5;%0.5;
L=1000;
%%%%%%%%%%%%%%%%%%%%%%%
% tfutr = 60e-3;
% w = 20e-3;
% x1=rectpuls(t-tfutr,w/2);
%%%%%%%%%%%%%%%%%%%%%%
% x1=square(t);
%%%%%%%%%%%%%%%%%%%%%
% x1= oscillator('Sawtooth',2,10); %2 second sawtooth at 440 Hz.
x1 = sin(2*pi*5*t);%+sin(2*pi*6*t)+sin(2*pi*8*t);
x=x1;%
% x=fftshift(x1);

fig=figure,plot(t,x');
xlabel('Time in seconds')
ylabel('Amplitude in volts')
title('Sinusoidal signal')
fig.Color=[1,1,1]
fig.Position=[1 1 1366 691]
% %%%%%%
% fig=figure;imshow(resIm,[]), title('original and ifrft on one plot')
% fig.Color=[1,1,1];
% fig.Position= [1 1 1366 691];
%%%%%

Faf = frft(x,a);
Fafs=abs(fftshift(Faf));
Mk=Fafs(1:length(Faf)/2);
fig=figure,subplot(1,2,2),plot(Mk);grid on,
xlabel('Frequency in Hz')
ylabel('Magnitude')
title('Transform using FrFT with alpha  = 2*pi')
fig.Color=[1,1,1]
fig.Position=[1 1 1366 691]
n = 2^nextpow2(L);
Y = fft(x,n);
P2 = (abs(Y)/n);
P1 = P2(:,1:n/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);  
f = (1:length(Y))*100/length(Y);
subplot(1,2,1);plot(0:(fs/n):(fs/2-f/n),P1(1,1:n/2)),hold on
xlabel('Frequency in Hz')
ylabel('Magnitude')
title('Transform using FFT ')

% % % % % % 
%  pp=[];
% fig=figure,
% xlabel('time period')
% ylabel('Amplitude in volts')
% title('obtained Sinusoidal signal using iFrFT')
% fig.Color=[1,1,1]
% fig.Position=[1 1 1366 691]
%  for k=1.01:0.005:1.03
    Faf1= (frft(Faf,-3*a/0.5));
    kk=ifft(Y);
    pk=real(Faf1);
%     plot((x1'),'b');hold on
%     plot(kk(1:length(t)),'k'),hold on
    figure,plot(pk,'r');
    %%%%%%%%%%%%%
%     cc=corrcoef(x1',pk);
      xlabel('time period')
 ylabel('Amplitude in volts')
 %title(['correlation=',num2str(cc(1,2))])

title('obtained Sinusoidal signal using iFrFT with alpha = -9pi ')
    %legend('inverse FFT','inverse FrFT - scaled')
%     legend('original','inverse FFT','inverse FrFT - scaled')
% %     pause(2)
%     hold off
 
  pp=[pp,cc(1,2)];
%  end
fig=figure,plot(t,pp)
title('correlation for input and obtained sinusoid')
fig.Color=[1,1,1]
fig.Position=[1 1 1366 691]

% % % % % % 


% pp=[];
% fig=figure,
% xlabel('time period')
% ylabel('Amplitude in volts')
% title('obtained Sinusoidal signal using iFrFT')
% fig.Color=[1,1,1]
% fig.Position=[1 1 1366 691]
% for k=1.01:0.005:1.30
%     Faf1= (frft(Faf,-k*a/1.00));
%     kk=ifft(Y);
%     pk=real(Faf1);
% %     plot((x1'),'b');hold on
%     plot(kk(1:length(t)),'k'),hold on
%     plot(pk,'r');
%     %%%%%%%%%%%%%
%     cc=corrcoef(x1',pk);
%     title(['correlation=',num2str(cc(1,2))])
%     
%     legend('inverse FFT','inverse FrFT - scaled')
% %     legend('original','inverse FFT','inverse FrFT - scaled')
%     pause(0.5)
%     hold off
%     pp=[pp,cc(1,2)];
% end
% fig=figure,plot(pp)
% title('correlation for input and obtained sinusoid')
% fig.Color=[1,1,1]
% fig.Position=[1 1 1366 691]

% % % % % % 
% pp=[];
% fig=figure,
% xlabel('time in seconds')
% ylabel('Amplitude in volts')
% title('obtained Sinusoidal signal using iFrFT')
% fig.Color=[1,1,1]
% fig.Position=[1 1 1366 691]
% %  k=1.01:0.005:1.30
%     Faf1= (frft(Faf,-1.01*a/1.00));
%     kk=ifft(Y);
%     pk=real(Faf1);
% %     plot((x1'),'b');hold on
%     plot(kk(1:length(t)),'k'),hold on
%     plot(pk,'r');
%     %%%%%%%%%%%%%
%     cc=corrcoef(x1',pk);
% %     title(['correlation=',num2str(cc(1,2))])
%     
%     legend('inverse FFT','inverse FrFT - scaled')
% %     legend('original','inverse FFT','inverse FrFT - scaled')
% %     pause(2)
%     hold off
%     pp=[pp,cc(1,2)];
% 
% fig=figure,plot(pp)
% title('correlation for input and obtained sinusoid')
% fig.Color=[1,1,1]
% fig.Position=[1 1 1366 691]
% % % % % % % % % 

% pp1=[];
% fig=figure,
% xlabel('Frequency in Hz')
% ylabel('Magnitude')
% title('obtained Transform using FFT and FrFT')
% fig.Color=[1,1,1]
% fig.Position=[1 1 1366 691]
% for k=1.01:0.005:1.30
% 
%     pk1= frft(pk,k*a);
%     pks=abs(fftshift(pk1));
% %     figure
%     plot(Mk)
%     hold on;
%     Mk1=pks(1:length(pk1)/2)
%     plot(Mk1);
%     cc1 = corrcoef(Mk,Mk1);
%     title(['correlation1=',num2str(cc1(1,2))])
%     legend('Frft of original ','Frft of retrieved ')
%     pause(0.5);
%     hold off;
%     pp1=[pp1,cc1(1,2)];
% end
% fig=figure, plot(pp1)
% title('correlation for transform and obtained transform')
% fig.Color=[1,1,1]
% fig.Position=[1 1 1366 691]
% 
% 
% 
% % figure,plot(real(Faf1))
% % figure,plot(imag(Faf1))
