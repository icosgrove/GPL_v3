function [quiet_whiten, quiet_fft, quiet_base, noise_floor, blocked, baseline0] = GPL_quiet(sp,sp_whiten,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function GPL_quiet will identify the signal absent, or quiet,
% portions of the window. The outputs are:
% quiet_whiten: Signal absent portions of the whitened spectrogram 
% quiet_fft: Signal absent portions of original spectrogram
% quiet_base: Mean noise level of entire slate
% noise_floor: Minimum noise value for the current window
% blocked: A vector where '0' is signal presence and '1' is signal absense,
% one value for each time bin.
% baseline0: Energy sums for each time bin, used for identifying signal
% presence/absence.

% Written by Tyler Helbe
% Tested, documented, and cleaned by Ian, but not fundamentally modifed
% 02/03/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Create signal present slate and baseline energy sum for each tiem bin.

% Create restricted energy sum frequency range, if necessary. Use whens
% umming energy values over signal only slate. 
low_off = parm.FreqBinLo - parm.FreqBinLo + 1;
high_off = parm.FreqBinHi - parm.FreqBinLo + 1;


% Take the RSS value for each time bin, project it to all frequency bins,
% and normalize whitened spectrogram over it.
u = sp_whiten./(ones(parm.NumFreqBins,1)*sum(sp_whiten.^2).^(1/2));
% Heble et al 2012 eq 7

% Take the RSS value for each frequency bin, project it to all time bins 
% and normalize the whitened spectrogram over it.
y = sp_whiten./(sum(sp_whiten'.^2).^(1/2)'*ones(1,parm.NumTimeBins));
% Heble et al 2012 eq 8

% The whitening step here was NOT in Helble et al. 2012
u = whiten_matrix(u')'; % Whiten normalized matrix over time
y = whiten_matrix(y); % Whiten normalized matrix over frequency

% Create the signal only slate using time/frequency exponents.
bas = abs(u).^parm.xp1.*abs(y).^parm.xp2; 
% Heble et al 2012 eq6

% Take the column-sum of the square of the signal only slate. Low_offf and
% high_off can restrict which rows are summed.
baseline0 = sum(bas(low_off:high_off,:).^2);






%% Find the quiet portions of the slate

% Create quiet slate starting with the whitened spectrogram from GPL_whiten
qs = sp_whiten; 

% At the end: correlate indices of quiet portions back to original
% spectrogram
quiet_fft = sp; 

lks = 0;

% Loop through time bins looking for singal absense. 
while lks < parm.NumTimeBins 
    
    % Recreate the signal present energy sum as done previously, but no
    % whitening. (Background noise not as heavily reduced as baseline0).
    u = qs./(ones(parm.NumFreqBins,1)*sum(qs.^2).^(1/2));
    y = qs./(sum(qs'.^2).^(1/2)'*ones(1,parm.NumTimeBins));
    bas = abs(u).^parm.xp1.*abs(y).^parm.xp2; 
    baseline1 = sum(bas.^2); 

    % Locate signal absent time bins (energy sums below noise ceiling).
    ks = find(baseline1 <= parm.noise_ceiling);
    
    % Set the quiet slate to be whitened spectrogram of only signal absetn
    % time bins. Same indices for quiet_fft but on original spectrogram.
    qs = qs(:,ks);
    quiet_fft = quiet_fft(:,ks);

    % Find number of signal absent time bins that were identifed
    lks = length(ks);  
 
    % Find number of times to extend random variable vector 
    xtnd = ceil(parm.NumTimeBins/length(ks));
    dex =[];
    
    % Create a vector of random integers between 1 and number of signal
    % absent time bins. Repeat that the number of times set by 'xtnd'
    for k = 1:xtnd 
        dex = [dex,GPL_shuffle(lks,lks)];                                     
    end

    % Restrict 'dex' back down to overall number of time bins
    dex = dex(1:parm.NumTimeBins); 
    
    % Mix up the quiet regions using randomness created by dex.
    qs = qs(:,dex);
    quiet_fft = quiet_fft(:,dex); 

end

% Tyler originally only saved 60 bins for quiet_fft, which due to the
% randomness will change the background noise level measurement (occurs
% later in the algorithm) with every run of the algorithm.
quiet_fft = quiet_fft(:,:);
%%% No signal presence removed %%%%


%% First Pass through Quiet times

% Whiten the baseline energy values
[b0,quiet_base] = whiten_vec(baseline0');
% b0: mean-reduced vector
% quiet_base: mean noise value for baseline0

% Normalize mean-reduced baseline vector by it's mean noise
baseline0 = b0'/quiet_base; 

% Identify the noise floor of the slate, the lowest energy value after all
% whitening, normalizing, and mean calculations. 
noise_floor = -min(baseline0); 


% Locate signal absent time bins (mod. energy sums BELOW the noise ceiling
% multiple). 
ks = find(baseline0 <= parm.noise_ceiling*noise_floor);

% Save 'qs' as signal absent time bins from the whitened spectrogram.
qs = sp_whiten(:,ks); 





%% Second Pass to idenify lower threshold calls in busy call environments

% Identify the number of signal absent columns found in first pass.
[~,tz2] = size(qs); 
% Find number of times to replicate the signal absent bins
rep = ceil(parm.NumTimeBins/tz2); 
qh = [];

% Duplicate the signal absent slate until it is longer than number of
% overall time bins.
for k = 1:rep 
    qh = [qh,qs];  
end

% Restrict down to overall number of time bins, 
qs = qh(:,1:parm.NumTimeBins);


% Repeat the RSS normalization, whitening, and exponent manipulation on the
% signal absent slate to extract signals not initally found.
u = qs./(ones(parm.NumFreqBins,1)*sum(qs.^2).^(1/2));
y = qs./(sum(qs'.^2).^(1/2)'*ones(1,parm.NumTimeBins));
u = whiten_matrix(u')';
y = whiten_matrix(y);
bas = abs(u).^parm.xp1.*abs(y).^parm.xp2; 

% Find column sum for second pass and whiten
baseline1 = sum(bas(low_off:high_off,:).^2); 
[b1] = whiten_vec(baseline1'); 


% Normalize the second pass baseline values by the mean noise of the FIRST
% pass. (or mean noise of the entire slate)
baseline1 = b1'/quiet_base; 

% Get rid of duplicate portions, now back to original size of signal absent
% slate.
baseline1 = baseline1(1:tz2);

% Locate baseline energy values ABOVE the noise ceiling (second pass signal
% presence).
ks_refine = find(baseline1 > parm.noise_ceiling*noise_floor);

% REMOVE the signal PRESENT bins that were just found.
ks = ks(setdiff((1:tz2),ks_refine));


%% Correlate signal absent indices back to the spectrograms

% Block out signal ABSENSE with a '1', presence with a '0'
blocked = zeros(parm.NumTimeBins,1); 
blocked(ks) = 1; 

% Locate signal absent time bins of the whitened spectrogram
qs = sp_whiten(:,ks); 

% Replicate the quiet portions of whitened sp so it can be as long as the
% original spectrogram
[~,tz2] = size(qs); 
rep = ceil(parm.NumTimeBins/tz2); 
qh = [];
for k = 1:rep
    qh = [qh,qs];
end

% Signal absent bins of whitened spectrogram (with duplication to fill it
% in).
quiet_whiten = qh(:,1:parm.NumTimeBins);
