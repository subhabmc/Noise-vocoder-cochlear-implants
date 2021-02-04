%% Clear workspace
close all
clear
clc

%% Initialize variables
f1     = 100;
f2     = 8000;
thresh = .4;
c      = 500;
Fs     = 10000;

% load the word captured
x1 = audioread('please.wav') ;
x=x1(:,1);
% Electrode center frequencies
CI_12 = logspace(log10(f1), log10(f2), 13);

% Random numbers with normal distribution
Noise = normrnd(0, 1, size(x));

% Matrix to store the filter channels
y_noise = zeros(length(Noise), length(CI_12)-1);
y_12    = zeros(length(x), length(CI_12)-1);
a       = zeros(3, length(CI_12)-1);
b       = zeros(3, length(CI_12)-1);

%% 1.1 Filtered noise

figure(1)
plot(Noise)
set(gca,'FontSize',12)
title('Noise signal')
print -depsc2 Plot1.eps
figure(2)
for i = 1:length(CI_12)-1
    hold on
    [b(:,i), a(:,i)] = butter(1, [CI_12(i) CI_12(i+1)]/Fs, 'bandpass');
    y_noise(:, i) = filter(b(:,i), a(:,i), Noise);
    y_12(:, i) = filter(b(:,i), a(:,i), x);
    plot(y_noise(:, i))
end
hold off
set(gca,'FontSize',12)
title('White Noise Filter Bank of each filter channel of a 12-electrode CI')
legend
print -depsc2 Plot2.eps

%% 1.2 Extract envelopes

% Implementation of Hilbert Transformation for the filtered speech signal
env = abs(hilbert(y_12));

figure(3)
plot(env)
title('Envelopes of the filtered speech signal')
set(gca,'FontSize',12)
legend
print -depsc2 Plot3.eps

%% 1.3 Add dynamic compression

% Compression
env_compressed = log10(1 + c*env)/log10(c + 1);

% Clipping
env_compressed( env_compressed > 1 ) = 1;
env_compressed( env_compressed < thresh ) = 0;

figure(4)
plot(env_compressed)
set(gca,'FontSize',12)
title('Dynamic Compression of Envelopes')
legend
print -depsc2 Plot4.eps

%% 1.4 Join the channels

% Modulation of band-limited noise 'y_noise' with the respective envelope 
% from the original speech signal 'env'
M = env_compressed .* y_noise;

% Summing all channels
Y = sum(M, 2);

% Listen the sum
sound(Y)
pause(2)

% Calculation of the single sided amplitude
SSAS_12 = Spectrum(Y);

% Single sided plot
figure(5)
stem(SSAS_12, "filled")
set(gca,'FontSize',12)
title('Single sided amplitude')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
grid on
print -depsc2 Plot5.eps
% Spectrogram of the joined channels
figure(6)
spectrogram(Y, round(Fs*10e-3), round(Fs*5e-3), Fs, 2*Fs, 'reassigned','yaxis')
title('Spectrogram of 12-electrode CI for word')
print -depsc2 Plot6.eps
