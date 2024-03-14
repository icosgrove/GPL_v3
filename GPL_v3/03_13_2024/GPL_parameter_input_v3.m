function [parm] = GPL_parameter_input_v3(PARAMS,xwav_struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_parameter_input_v3 is a script to set fixed input parameters for the
% GPL detector algorithm. Some parameters are deterined by xwav header 
% information, others are manual input Some of these parameters will change 
% during the GPL process. 

% Written by Ian Cosgrove 01/26/2024
% Lsst Updated 01/26/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Sample Frequency [samples/s]: Value straight from xwav header
parm.SampleFreq = PARAMS.ltsahd.sample_rate(1);

% Number of samples per xwav 
%%% nxwav works for HARP data, doesnt work for
% xwavs saved manually from Triton, is there a way to distinguish? %%%
parm.nrec = xwav_struct.xwav1.TotalSamples/xwav_struct.xwav1.nxwav;
% 


% Helble et al 2012 eq 6 (nu1 and nu2)
% GPL Algorithm Parameters
parm.xp1 = 2.5; % Exponent for u: Normed over frequency
parm.xp2 = 0.5; % Exponent for y: Normed over time
parm.whiten = 1; % 1: Whiten in time only. 2: Whiten in time and frequency
parm.white_x = 1.5; % Enhanced mean noise scale factor
parm.NumLoops = 5; % Max number of loops each window is processed for calls



% Desired Frequency Range to be Processed [Hz]
parm.freq_lo = 30;
parm.freq_hi = 90;

% Call duration Minimum/Maximum [s]
parm.MinCallDuration = 0.5;
parm.MaxCallDuration = 10;

    
    % Spectral Parameters
    parm.fftl = 2000; % FFT Length
    parm.fftOverlapPercent = 90; % FFT overlap percentage
    parm.fftOverlap = parm.fftl - (parm.fftOverlapPercent*1e-2)*parm.fftl; % FFT bin step size



    % Number of Time bins per window (Number of columns in spectrogram)
    parm.NumTimeBins = floor((parm.nrec - parm.fftl)/parm.fftOverlap) + 1;
    





% Create Frequency Range in bins based on FFTL and sampling rate. 
freq_range = 0:(parm.fftl/2)/parm.fftl*parm.SampleFreq; % Define frequency range from FFTL

% Identify Frequency bin range
[~,parm.FreqBinLo] = min(abs(freq_range - parm.freq_lo));
[~,parm.FreqBinHi] = min(abs(freq_range - parm.freq_hi));
parm.NumFreqBins = parm.FreqBinHi - parm.FreqBinLo + 1;

% Frequency bin range for summation (Currently unchanged from original
% range)
parm.SumFreqBinLo = parm.FreqBinLo;
parm.SumFreqBinHi = parm.FreqBinHi;



% Detection/Contour Threshold Parameters
% One-Dimensional energy summation thresholds
parm.noise_ceiling = 50; % Upper limit for backgorund noise, multiple of noise floor
parm.thresh = 200; % Call threhsold: minimum multiple of noise floor for a call to pass

% Two-Dimensional Energy Threshold for each pixel in the spectrogram
parm.ContourCutoff = 6; % Minimum value for time,frequency coordinate to be included to the call contour
parm.MinIslandSize = 10; % Minimum number of pixels above cutoff for a detection to be considered a contour





%% Measurement On/Off Switches
parm.waveform = 1; % Extract the original time domain waveform of the call
parm.template = 1; % Extract call contour and measurements (0: only start/end time is recorded)
parm.cm_on = 1; % Extarct the call contour 
parm.cm_max_on = 1; % Extract the strongesr contour of the call
parm.cm_max2_on = 1; % Extract the second strongest contour of the call
parm.measurements = 1; % Obtain measurements of the contours
parm.filter = 0; % Run primary filters %% EXTRA parm to apply filter?


    % Measurement on/off swtiches
    parm.measure.cm_freq_limits = 1; % Measure frequency start/end for Cm
    parm.measure.cm_max_freq_limits = 1; % " " for cm_max
    parm.measure.cm_max2_freq_limits = 1; % " " for cm_max2

    parm.measure.cm_freq_bandwidth = 1; % Measure freq bandwidth for cm
    parm.measure.cm_max_freq_bandwidth = 1; % " " for cm_max
    parm.measure.cm_max2_freq_bandwidth = 1; % " " for cm_max2

    parm.measure.cm_duration = 1; % Measure call duration for cm
    parm.measure.cm_max_duration = 1; % " " for cm_max
    parm.measure.cm_max2_duration = 1; % " " for cm_max2

    parm.measure.cm_slope = 1; % Measure slope for cm
    parm.measure.cm_max_slope = 1; % " " for cm_max
    parm.measure.cm_max2_slope = 1; % " " for cm_max2

    parm.measure.spec_noise = 1; % Measure background noise of the call's xwav
    parm.measure.spec_rl = 1; % Measure the calls recieved level relative to background noise




%% Preliminary Filtration Switches and Values:

% On/off: Merge calls within time bin radius of another call (handle
% intra-call energy drops).
parm.filter_parm.switches.AdjacentCallMerger = 1;

% Set bin number for merging:
parm.filter_parm.values.AdjacentCallBinNum = 3;



    parm.filter_parm.values.cm_max_duration = 2; % Minimum time length (time bins) for cm_max
    parm.filter_parm.values.cm_max2_duration = 15; % " " for cm_max2

    parm.filter_parm.values.cn_max_bandwidth = 3; % Minimum frequency bandwidth (frequency bins) for cm_max
    parm.filter_parm.values.cn_max2_bandwidth = 3; % " " cm_max2

    parm.filter_parm.values.cm_max_slope_lower = -1.8; % Minimum slope for cm_max
    parm.filter_parm.values.cm_max2_slope_lower = -1.8; % " " cm_max2
    parm.filter_parm.values.cm_max_slope_upper = 0; % Maximum slope for cm_max
    parm.filter_parm.values.cm_max2_slope_upper = 0; % " " cm_max2


% Plot the detection window (will stop the detection process after each
% window containing a call until user manually progresses.
parm.plot = 0;


% Island radius: Number of bins to check around islands for others (connect
% the dots for contours)
parm.island_radius = 1; %% CHANGE with better contour connection emthods



%% Save Output

save('sei_whale_90_30_parm.mat','parm');

end % Function

