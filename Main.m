%% Prepare working space
clear
close all
clc

fa=100;     %Apical electrode frequency
fb=8e3;     %Basal electrode frequency

%Number of electrodes
n1=12;      

%Logspace for electrode centre frequencies
al=log(fa)/log(10);
bl=log(fb)/log(10);
f1=logspace(al,bl,n1+1);


fs=20e3;  %Sampling frequency
c=500;

%Initializing
% 3-Electrode CI
a1=zeros(length(f1)-1,3);
b1=zeros(length(f1)-1,3);
yfs1=0;

[y,Fs] = audioread('please.wav') ;
noise = 1*randn(1,length(y));
Env1=zeros(length(f1)-1,length(y));
Env1c=zeros(length(f1)-1,length(y));
Env1m=zeros(length(f1)-1,length(y));
Env1mr=zeros(1,length(y));
%Filter bank implementation, frequency response, output, spectra and spectrograms
for i=1:length(f1)-1
order    = 1;
fcutlow  = f1(i);
fcuthigh = f1(i+1);
x=[fcutlow,fcuthigh]/(fs/2);
%Butterworth filter
[b1(i,:),a1(i,:)] = butter(order,x, 'bandpass');
[h,f] = freqz(b1(i,:),a1(i,:),2*fs,fs);
%Input audio file

y=y(:,1);
%Filtered audio file (Individual filters)
yf = filter(b1(i,:),a1(i,:),y);
Env1(i,:)=abs(hilbert(yf));
Env1c(i,:)=log10(1+c*Env1(i,:))/log10(1+c);
Env1c(Env1c > 1) = 1;
Env1c(Env1c < 0.1) = 0;
Env1m(i,:)=Env1c(i,:).*noise(1,:);

%Resultant from filters
yfs1=yfs1+yf;
%Butterworth filter frequency response
figure (i);
plot(1:length(yf),Env1m(i,:));

Env1mr(1,:)=Env1mr(1,:)+Env1m(i,:);
end
sound(Env1mr,Fs);

ffsig3=fftsignal(Env1mr);
L=length(Env1mr);
fsig = Fs*(0:(L/2))/L;
figure (7);
stem(fsig/1000,ffsig3);
title('Single-sided amplitude spectrum of 12-Electrode CI');
xlabel('Frequency (kHz)');
ylabel('Amplitude');
print -depsc2 Plot1.eps

%12-Electrode CI
figure (11);
colormap(jet);
spectrogram(Env1mr,round(Fs*10e-3),round(Fs*5e-3),Fs,2*Fs)
title('Spectrogram of 12-Electrode CI for given input word');
ax = gca;
ax.YDir = 'normal';
print -depsc2 Plot12.eps


