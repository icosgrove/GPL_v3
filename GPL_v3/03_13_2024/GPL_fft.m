function [sp] = GPL_fft(data,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_fft take the Fast Fourier Transform of the current window of audio
% data with a given FFTL and overlap set by the user in the input
% parameters and create a spectrogram.
% Modified by Ian Cosgrove 02/02/2024

% Inputs: 
% data: Current window of audio data
% parm: Input parameters
% Outputs:
% sp: Spectrogram of the current window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Allocate and create Hammming WIndow
win = hamming(parm.fftl); % Hamming Window
sp = zeros(parm.NumFreqBins , parm.NumTimeBins); % Spectrogram Window

x = data;
% Transpose to column vector if necessary
[x1,x2]=size(x);    
if x2>x1
    x=x'; 
end 

%% Take FFT and Create Spectrogram
for j = 1:parm.NumTimeBins % Loop through all time bins
    
    % Start/End Samples 
    start = (j-1)*parm.fftOverlap + 1; % Moving start point for FFT
    finish = start + parm.fftl - 1; % Moving end point for FFT
    
    % FFT and Spectrogram
    q = fft(x(start:finish).*win); % FFT
    sp(:,j) = abs(q(parm.FreqBinLo:parm.FreqBinHi)); % Write FFT coefficients to spectrogram

end

