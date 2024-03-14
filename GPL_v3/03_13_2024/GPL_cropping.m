function [base_out,newpad] = GPL_cropping(base_in,noise_ceiling,thresh)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function GPL_cropping will identify detections for the current window
% using, time bin column energy sums, the noise celling parameter, and
% detection threshold parameter. The function will search for the highest
% energy portions of the window, stradled by signal absense, and save them.
% It will then remove that detection, and search for the next highest
% energy valued section, repeating until all units above the detection
% threshold are found.

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian, but not fundamentally modified
% 02/03/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Allocate outputs
iflag = 0;
newpad = [];
base_out = base_in;

% Loop until all units above threshold have been identified
while iflag == 0
    
    % Check if the maximum energy value is less than detection threshold,
    % and if so end unit detection.
    [k1,k2] = max(base_out); 
    if k1 < thresh 
        return 
    end
    
    % Locate signal absent time bins (Below noise ceiling)
    js = find(base_out < noise_ceiling);  
    
    
    % Left straddle:
    % Locate signal absent bins that occur BEFORE the maximum energy bin (in time)
    j1 = find(js < k2);
    
    % Locate the index closest signal absent bin BEFORE the highest energy bin
    if isempty(j1) == 1 
        js1 = 1;
    else
        js1 = js(j1(end)); 
    end
    
    
    % Right straddle:
    % Locate signal absent bins that occur AFTER the maximum energy bin (in time)
    j2 = find(js > k2); 
    
    % Locate the index of the closest signal absent bin AFTER the maximum energy bin
    if isempty(j2) == 1 
        js2 = length(base_out);
    else
        js2 = js(j2(1)); 
    end
    
    
    % Those indices found straddle the highest energy detection, or
    % indicate it's start and end. Save that as the detection.
    newpad = [newpad',[js1,js2]']'; 
    
    % Remove the highest energy detection, search for the signal absent
    % bins that surround the next highest one.
    base_out(js1:js2) = 0; 
    
    % Repeat until all units above threshold are identified.
end
