function [GPL_struct] = GPL_template(GPL_struct,sp,sp_whiten,start,finish,x,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_template will take each detection and determine where in time and
% frequency the energy spikes occured that originally flagged it as a call 
% during summation. This will idenitfy the contour of the call which is 
% used for later filtering. This function will also compress the contour 
% data and load it into GPL_struct. The original waveform of the call is 
% also computed here.

% Written by Tyler Helble 
% Tested, documented, and cleaned by Ian, but not fundamentally modified
% 02/10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m1,~] = size(start);

% Loop over each detection
for k = 1:m1 
    
    % Create whitened spectrogram with selected calls only
    base = sp_whiten(:,start(k):finish(k)); 

    
    % Create the call contours
    [cm, cm_max, cm_max2] = GPL_contour(base,parm.ContourCutoff,parm);

    
    
    if(parm.cm_on == 1)
    GPL_struct = GPL_sparse(cm,'cm',k,GPL_struct);
    end

    % This produces index values for what part of the call is cm, the
    % spectrogram values normalized an2 scaled to 2^16, the peak energy of cm,
    % and the size of cm in bins. All of this in a structure within GPL_struct

    if(parm.cm_max_on==1)
    GPL_struct = GPL_sparse(cm_max,'cm_max',k,GPL_struct);
    end
    if(parm.cm_max2_on==1)
    GPL_struct = GPL_sparse(cm_max2,'cm_max2',k,GPL_struct);
    end

    % The same process done with cm is repeated for cm_max and cm_max2

    %%% reconstruct with:
    %%%  cm = GPL_full('cm',k,GPL_struct);

    if parm.waveform==1 % set to 0 in blue whale D call parm file
     xv=x(1+(start(k)-1)*parm.skip:(finish(k)-1)*parm.skip);
     % Find the sample indices for the start/end of detection, and save the
     % samples of each detection to 'xv'


     [vec]=ww365(xv,cm_max,parm,sp(:,start(k):finish(k)));
     GPL_struct(k).cm_max_waveform=vec;
    end




end