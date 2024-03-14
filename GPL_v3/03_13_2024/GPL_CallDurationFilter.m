function [start,finish,times] = GPL_CallDurationFilter(mcalls,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_CallDurationFilter is adapated from Tyler Helble's GPL_prune, which
% removed detections that fall above or below the call duration
% parameters. A conversion is done to determine how many time bins
% correspond to seconds, and then calls with too little or too many time
% bins are removed. 

% Originally written by Tyler Helble as 'GPL_prune'
% Tested, documented, and cleaned by Ian, but not fundamentally modified
% 02/10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate the number of seconds per time bin
dt = parm.fftOverlap/parm.SampleFreq; 


%%% Potentially to be removed.
% Merge adjacent calls within allowed gap size  

[~,k2] = sort(mcalls(:,1)); % sort start indices and make index vector
start = mcalls(k2,1);
finish = mcalls(k2,2); 

k = find(start(2:end)-finish(1:end-1) < parm.filter_parm.values.AdjacentCallBinNum); 
% find calls within 5 bins of each other (overlapped), none in this loop
omit = length(k);
if omit > 0 % Combine adjacent calls 
  n=1:length(start); 
  start = start(setdiff(n,k+1));
  finish = finish(setdiff(n,k));
end
%%%


% Filter out calls shorter than minimum call length cutoff
cutoff_short = ceil(parm.MinCallDuration/dt) - 2; 
k = find(finish - start > cutoff_short); 
start = start(k);
finish = finish(k);

% Filter out calls longer than maximum call length cutoff
cutoff_long = ceil(parm.MaxCallDuration/dt)-2; 
k = find(finish - start < cutoff_long); 
start = start(k);
finish = finish(k); 


% Convert from bin number relative to window start to sample number
% relative to window start. Save final times after duration filtering
times = [(start-1)*parm.fftOverlap + 1,(finish)*parm.fftOverlap];
 
end
