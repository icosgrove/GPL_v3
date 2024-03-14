function sp_whiten = GPL_whiten(sp,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function GPL_whiten will create a whitened version of the original
% spectrogram. This is done by whitening the spectrogram over time by
% finding the frequency bin mean noise levels, reducing the spectrogram by
% those means, and normalizing it over the means. This may optionally be
% done over frequency by doing the same process with individiual time
% bins. The spectrogram then undergoes a 2-D convolution which has the
% affect of blurring the spectrogram minimally.

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian, but not fundamentally changed.
% 02/02/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Find the frequency mean noise and mean-reduced spectrogram
[sp_whiten,~,mu] = whiten_matrix(sp,parm.white_x);

% mu: Vector of mean noise (signal-absent) level for each frequency bin
% sp_whiten: Spectrogram with mean noise across each frequency bin
% subtracted from original values. 




%% Create whitened spectrogram (whitened in time)

% Normalize the mean-reduced spectrogram by its mean noise values (One mean
% for each frequency bin expanded to cover all time bins).
sp_whiten = abs((sp_whiten./(mu*ones(1,parm.NumTimeBins))).')';


%% Optional: Whiten in frequency (on top of time)
if parm.whiten == 2 
    
    % Take the frequency-whitened spectrogram and produce mean noise levels
    % for each of it's time bins, and use that to reduce it further.
    [sp_whiten,~,mu] = whiten_matrix(s-p_whiten');
    
    % Normalize the spectrogram over it's time bin mean noise values. (One
    % mean for each time bin expanded to cover all frequency bins).
    sp_whiten = abs(sp_whiten'./(ones(parm.nfreq,1)*mu'));
 
end



%% Perform 2-D convolution on the whitened spectrogram

% Create kernel matrix
cross = ones(3,3);
cross(2,2) = 4;
cross(1,3) = 0;
cross(3,1) = 0;
cross(1,1) = 0;
cross(3,3) = 0;
cross = cross/8;

% Perform the convolution (smoothing)
sp_whiten = conv2(sp_whiten,cross,'same');
