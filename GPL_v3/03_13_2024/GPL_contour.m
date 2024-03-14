function [cm,cm_max,cm_max2] = GPL_contour(bas0,cutoff,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_contour will take each detection individually, perform manipulations
% and identify where in time and frequency it's energy spikes occurred. The
% contour 'cm' is all of these high-energy pixels. Pixels are separted from
% each other by regions of low energy, creating islands. 'cm_max' is the
% island with the most energy, 'cm_max2' is the island with the second most
% amount of energy. Measurements are later taken of cm_max only, in order
% to focus on the primary contour of the call which is likely the signal of
% interest itself without anything extraneous. 

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian but not fundamentally modified.
% 02/10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Allocate 
% Add a one bin buffer around the call
[sz1,sz2] = size(bas0); 
bas = zeros(sz1+2,sz2+2); 
bas(2:end - 1,2:end - 1) = bas0; 

% Allocate Contour matrices 
cm = 0*bas;
cm_max = cm;
cm_max2 = cm;
ks0 = 0;




%% Apply whitening and normalization to the call sp data
% Reshape the call spectrogram into a single vector 
[sz1,sz2] = size(bas); 
basgot = reshape(bas,sz1*sz2,1); 

% Whiten the vector and normalize by it's mean noise level
[~,mu] = whiten_vec(basgot); 
qre = bas/mu - 1; 

% If not provided in input args, define the pizel energy cutoff
if nargin == 1 
    cutoff = 10; 
end




%% Identify Islands

% Define minimum number of pixels for a detection to be considered a contour.
a_min = parm.MinIslandSize;

% Locate indices of whitened/normalized pixel levels above the energy
% cutoff. These are pixels considered as part of the contour
k = find(qre > cutoff); 


% Identify islands of the contour if it passes minimum length filter
if length(k) > a_min 
    
    % Set passed pixel indices to '1' in a window the size of the call
    msk0 = zeros(sz1,sz2); 
    msk0(k) = 1; 
    
    % Perform a 2-D convolution using the central part on the passed pixels
    dm = msk0; 
    dm = conv2(dm,ones(3,3)/9,'same');
    
    % get rid of low-lying pixels after convolution
    dm(dm > 2/9) = 1; 
    dm(dm < 1) = 0; 
    
    % Apply the convolution (with some pixels on the edges removed in
    % previous step) to the original pixels. This has the effect of adding
    % some extra pixels above cutoff around the call.
    dm = dm + msk0; 

    % Locate new set of pixel indices
    kk = find(dm); 

    
    %%%%%%%%%%%%%%%%%% UPDATED VERSION ENDS HERE - Ian 03/13/2024
    %% Work in progress
    
    
    % Locate islands in the call. Choice of expanded radius available
    if parm.island_radius == 1
        chain1 = island1x(kk,sz1);
        %chain2 = island1x_OLD(kk,sz1,sz2);
    else
        chain = radial_island(kk,sz1,sz2,parm);
    end
    % Output of island1x: 'chain' containing the length of each island
    % followed by the indices of the island over the individual call window
 

    
    [area,t,col] = island21(bas,chain1); 

    
    
 if max(area)>a_min % check if at least one island is greater than 10 bins. Max: 387
 score=1; %%%% at least one viable contour exists, proceed
   ks=find(area>0.05*max(area)); % find islands with size > 19.35 bins only 1 & 4 

 msk=zeros(sz1,sz2);  % allocate 59x75
 [dummy,ks0]=sort(t); % ks0 is the indices of the sorted island values, smllest to largest. 
 % Island 6 is the smallest and island 1 is the largest
     msk(chain(col(ks0(end))+1:col(ks0(end)+1)-1))=1;
     % First find ks0(end) or the index of the largest call (1)
     % col(1) will find starting index of the largest call
     % Add one to that to account for the separation element in chain
     %   - Each island's bins are separated by the length of the next
     %   island
     % col(ks0(end)+1)-1) will find the end index of the largest island
     % chain(" ") will pull the bins for the island
     % msk(" ")=1 will set the largest island to 1s. 
     
     
     
     
     cm_max=msk.*bas;  %%%% Single contour of biggest area in the slate
     % multiplying will set all the bins that are not part of the largest
     % island to 0. 
     
     % cm_max/largest contour/largest feature matrix is created, which is
     % the strongest energy island of the spectrogram that are above the
     % noise threshold. For this case the largest island in # of bins
     % happened to be the strongest contour, but th real determining factor
     % of the contour is it's energy intensity and not raw size.
     
     
     if(length(ks0)>1) % if there is more than 1 contour larger than 10 bins, proceed 
      msk=zeros(sz1,sz2); % allocate 59x75
        msk(chain(col(ks0(end-1))+1:col(ks0(end-1)+1)-1))=1; % find second strongest island
        cm_max2=msk.*bas;  %%%% Single contour of biggest area in the slate
        % write second strongest contour to cm_max2
     end
     
     % The same process si repeated in finding the bins for the eccond
     % strongest contour, if it exists 
 
     msk=zeros(sz1,sz2); % sllocate 59x75
     for dummy=1:length(ks) % loop over the contours above 5% cutoff only. Only 2 islands for our example
     msk(chain(col(ks(dummy))+1:col(ks(dummy)+1)-1))=1; % combine the contours together 
     end
     cm=msk.*bas; % write all of the strongest contours to 'cm', could be more than 2
 end
end

cm=cm(2:end-1,2:end-1); % get rid of extra bin, size is back to 57x73
cm_max=cm_max(2:end-1,2:end-1); % " "

if(length(ks0)>1)
cm_max2=cm_max2(2:end-1,2:end-1); % if the second strongest contour exists, reduce back to 57x73
end
a=1;




