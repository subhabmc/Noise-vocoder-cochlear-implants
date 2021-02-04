% 2. Calculation of the Single Sided Amplitude Spectrum (positive half) of 
% the Signal using Fast Fourier Transformation

% The spectrum function takes the existed signal and returns the Single
% Sided Amplitude Spectrum
function [SSAS] = Spectrum(signal)

% cutting off negatives with absolute and adding up the positives
SSAS = 2*abs(fft(signal)/length(signal));
SSAS = SSAS(1:length(signal)/2+1);
% keep offset and last value unaffected
SSAS(1) = SSAS(1)/2;                             
SSAS(end) = SSAS(end)/2;


end