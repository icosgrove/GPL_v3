function [ks] = base3x(baseline)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function base3x.m will process vectors corresponding to a frequency
% bin over all time or a time bin over all frequency to find the indices
% where it is not at or near a maximum/minimum in energy.

% The function sorts the energy values smallest to largest, splits that in
% half and subtracts the smaller energy from the higher energy values. The
% smallest difference found is then attributed to some index in the bottom
% and top half of the values. The energy values above and below the cutoff
% are removed, leaving only the 'middle' values. Those middle indices are
% attributed back to there original location in the input vector and saved
% as the output to the function.

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian, but not fundamentally modified.
% 02/02/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Sort and Reshape the values

% Sort energy for each time bin over the current frequency bins into column vector
bs = sort(baseline)';

bs1 = length(bs); % Number of time bins
z = ones(bs1,1); 

% Make sure the number of time bins is even
if mod(bs1,2) == 1 
    bs = bs(1:end-1);
    bs1 = bs1-1;
end

% Split the sorted energy values at the middle and separate the two pieces
bs = reshape(bs,bs1/2,2);




%% Identify the middle energy value indices of original bins

% Identify the column containing the smallest difference between the
% separated sorted energy values
[~,b2] = min(diff(bs'));

% Identify the energy values associated with the smallest energy difference
% These indicate low and high energy cutoffs, the indices between which are
% of interest. 
low = bs(b2,1);
high = bs(b2,2);

% Remove indices below the energy minimum cutoff and above the maximum ctuoff.
% These are the original indices, not sorted.
z(baseline < low) = 0;
z(baseline > high) = 0;

% Set the output to be the original (non-sorted) indices containing energy
% values between the high and low cutoffs. 
ks = find(z);
