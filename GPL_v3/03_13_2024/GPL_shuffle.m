function [cr] = GPL_shuffle(nc,deck)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function GPL_shuffle will randomly shuffle the indices of the signal
% absent time bins.

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian, but not fundamentally modified.
% 02/03/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a vector of Gaussian random variables for number of signal absent
% time bins. Sort them from smallest to largest, and save original indices.
[~,cr] = sort(rand(1,deck));

% Save only indices specified by 'nc'. THe only use of this function has
% the number of signal absent time bins as both input arguments, so all the
% values will be saved. 
cr = cr(1:nc);
