function [area,t,col] = island21(bas,chain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Island21 will compute energy sums for each individual island. This is
% used to determine which island is the strongest and second strongest
% contour, which are the primary modes for filtration later

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian, but not fundamentally modified
% 02/14/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Code restructuring currently paused here - Ian 03/13/2024




area(1) = chain(1);
loc = 1; 
la = length(chain); % total number of indices for all islands
     
    for i=2:99999 % while idk
      loc=loc+area(i-1)+1; % Loop 1: loc=1+387+1 = 389
      
       if loc == la+1 % break loop is island index goes over the maxium index for all islands
         break
       end
       
      area(i)=chain(loc); % sets eecond index of area to qc(389). 
      % That corresponds to the length of the next island 
    end
% We end up getting an array containing the lengths of all the islands
% (number of indices of all the islands)
% area: [387 8 17 118 5 2]

  col(1)=1; % allocate column as just a 1
  
    for i=1:length(area) % loop over all islands
     col(i+1)=col(i)+area(i)+1; % this finds all the starting indices of the islands  
    end
    % col has the starting indices of the islands, diff(col) is the length
    % of the islands+1

    for i=1:length(area) % loop over all islands
     t(i)=sum(bas(chain(col(i)+1:col(i+1)-1))); 
    end
    a=1;
    % col(i)+1:col(i+1)-1 will select the range of indices of each island,
    % but have a bin padding on each side. So island 1 is index 1-389, and
    % this restriction gives us 2-388. 
    
    % qc(last step) will pull indices of this island
    % q(previous) wlll pull the whitened spectrogram values of the island
    % those values are summed and written to 't'
    
    % Outputs:
    % area: The lengths of each island
    % t: the summed spectrogram values of each island
    % col: the starting indices of each island

