function [chain] = island1x(orig_index,sz1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% island1x will locate islands of high energy pixels based on if a pixel
% lies directly next to (not including diagonals) another pixel. The final
% output is a vector containing the length of an island followed by it's
% indices.

% Written by Tyler Helble 
% Tested, documented, and cleaned by Ian, but not fundamentally modified.
% 02/10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocate
index_copy = orig_index; 
chain = [];
q = length(index_copy); 

% Loop over each individual pixel index
for q1 = 1:q 
    
    % Check if the pixel is already accounted for (NaN), this will
    % effectively skip to the next island once one is identified
    if isfinite(index_copy(q1)) == 1 

        % Load in current index to be checked (not necessarily consecutive
        % because indices can be acccounted for and NaN'd out)
        kk = index_copy(q1);
        k = []; 
        
        % Loop until there are no indices with high energy directly
        % adjacent to all the others in a built island
        while isempty(kk) == 0 
                
                % Append confirmed adjacent indices onto each other
                k = [k,kk'];
                
                % Clear from previous iteration
                surrounding_indices = [];  
                
                % Locate indices around the island, growing from just one
                % index, and eventually ch  ecking the idncies around the
                % entire island once the while loop is near termination.
                for i = 1:length(kk) 
                    
                    % Replace index being checked with NaN so it isn't
                    % re-checked
                    index_copy(index_copy == kk(i)) = nan; 
                    
                    % Create a vector of the indices directly above, below,
                    % to the right and left of the index being checked.
                    surrounding_indices = [surrounding_indices,[kk(i) + 1,kk(i) - 1,kk(i) + sz1,kk(i) - sz1]];

                end % End surrounding index loop
                
        
            % Save all immediate indices just found (grows as the island grows)
            kk = surrounding_indices; 
            
            % Identify if the surrounding indices are part of the same
            % island by checking for overlap with the original contour
            % indices.
            kk = intersect(kk,index_copy);
            

        end % End while loop: A full island has been identified
    
    % Append the length of the island and the island's indices to the
    % output
    chain = [chain,length(k),k];
    
    end % End index finite check
    
end % End loop over all indices


end


