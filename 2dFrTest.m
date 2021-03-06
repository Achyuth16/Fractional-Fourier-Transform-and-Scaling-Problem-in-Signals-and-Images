% function 2dFrTest

clc
clear all
close all
%%%%%%%%%%%%%%%%%
fs=255;
t=0:1/fs:1;
%0.5;
% x= sin(2*pi*rand(1)*t);
% x1 = sin(2*pi*2*t).^1+sin(2*pi*4*t)+sin(2*pi*6*t)+sin(2*pi*8*t);
%%%%%%%%%%%%%%%%%
Im = zeros(256);

for pp=1:30:255 
    Im(pp,:)= 128+128*sin(2*pi*rand(1)*t);
%     figure,plot(Im(pp,:))
end

fig=figure,imshow(Im,[]), title('original image')
fig.Color=[1,1,1]
fig.Position=[1 1 1366 691]

bb=fft2(Im);
fig=figure,subplot(1,2,1),imshow(log(abs(fftshift(bb))+1),[]);title('fourier transform of an image')
fig.Color=[1,1,1]
fig.Position=[1 1 1366 691]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(Im);
tmp=zeros(size(Im));
Phi=10;
for K1=1:m
        tmp(K1,:)= frft(Im(K1,:),Phi);      
end
FrM=zeros(size(Im));
for K2=1:n
    FrM(:,K2)= frft(Im(:,K2),Phi);   
end
subplot(1,2,2),imshow(abs(FrM),[]), title('frft of an image with alpha = 3pi/4')
fig.Color=[1,1,1]
fig.Position=[1 1 1366 691]

r=FrM;
g=Im;
bCh=zeros(256,256);
resIm(:,:,1)=r;
resIm(:,:,2)=g;
resIm(:,:,3)=bCh;

fig=figure,imshow(resIm,[]), title('pattern changes with transform angle')
fig.Color=[1,1,1]
fig.Position=[1 1 1366 691]

% for pp=1:30:255
%     figure,plot(Im(pp,:),'b'), hold on,plot(FrM(pp,:),'r')
% end    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   for K1=1:m
%        for k=1.01:0.005:1.30
        tmp(K1,:)= frft(Im(K1,:),-1.01*Phi);      
%        end
    end
FrM=zeros(size(Im));
for K2=1:n
    %for k=1.01:0.005:1.30
    iFrM(:,K2)= frft(tmp(:,K2),-3*Phi);   
    %end
end
fig=figure,imshow(real(iFrM),[]), title('ifrft of an image with alpha = -30')
fig.Color=[1,1,1]
fig.Position=[1 1 1366 691]

rCh=iFrM;
gCh=Im;
bCh=zeros(256,256);
resIm(:,:,1)=rCh;
resIm(:,:,2)=gCh;
resIm(:,:,3)=bCh;
fig=figure,imshow(resIm,[]), title('pattern chnages with change in angle')
fig.Color=[1,1,1]
fig.Position=[1 1 1366 691]

 
