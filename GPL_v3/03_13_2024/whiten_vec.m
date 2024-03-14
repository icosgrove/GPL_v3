function [spc,mu]=whiten_vec(sp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function whiten_vec will take an input vector, use base3x to identify
% the median noise indices, use that to find mean noise value of the
% vector, and finally create a mean-reduced version of original vector.

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian, but not fundamentally modified
% 02/03/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Identify median noise indices
[ks] = base3x(sp);

% Locate original noise values pertaining to median values, find mean noise
qs = sp(ks);  
mu = mean(qs,1);

% Create a vector with each entry reduced by the vector's mean noise.
spc = sp - mu;