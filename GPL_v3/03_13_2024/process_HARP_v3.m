%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process_HARP_v3 is the main script that will run the Generalized Power
% Law (GPL) Detector Algorithm on a set of acoustic data. The algorithm
% will parse through provided xwavs and perform statistical cmputations to
% extract cetacean vocalizations, or other non-biological signals of
% interest, from the data.

% This version was written by Ian Cosgrove with Joshua Jones starting
% 01/12/2024. The original versions (GPL_v2 (2020), GPL 2014) were created
% by Tyler Helble. See his publication: "A generalized power-law detection 
% algorithm for humpback whale vocalizations" (2012)
% (https://www.cetus.ucsd.edu/publications.html) for full documentation on
% the theory behind the GPL detector. 

% This script will perform the functions of loading the xwavs that are to
% be processed by the detector. The xwav headers will be extracted and the
% total amount of time will be separated into samller windows, these 
% windows will be processed one at a time in the script GPL_v3. After
% GPL_v3 and it's subfunctions have completed the detection process on a
% window, the output data will return to this script as GPL_struct, and it
% will be saved to the main output file. Once all windows have been
% processed, the final detection output file will be saved from here. 

% Documentation written by Ian Cosgrove
% Contact: icosgrove@ucsd.edu
% Last Updated: 03/13/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load parameter file and xwavs - Tyler Helble

% Load xwav files 
fprintf('\nSelect folder containing .x.wav to process\n');
    [filename, pathname] = uigetfile('*.x.wav'); % Select File
    cwd = pwd; 
    cd(pathname) % Set current directory to path containing xwav file
    file_dir = pwd; 
    addpath(pwd); 
    files = dir('*.x.wav'); % Variable files contains xwav data
    cd(cwd); % Set current directory back to current working directory


% Allocate array for detections
calls = [];

% End Helble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract xwav headers - Ian Cosgrove

% Pre-allocate
field_name = cell(length(files),1);

% Loop through each xwav in 'files' and extract xwav headers and audio info
for n = 1:length(files) 
    xwav_name = files(n).name; 
    PARAMS = getxwavheaders(file_dir,xwav_name); % Retrieve xwav data
    field_name{n} = sprintf('xwav%i',n);
    xwav_struct.(field_name{n}).julian_start_time = PARAMS.ltsahd.dnumStart; % xwav start time in julian time
    xwav_struct.(field_name{n}).year = PARAMS.ltsahd.year; % date string for start times
    xwav_struct.(field_name{n}).month = PARAMS.ltsahd.month; % " "
    xwav_struct.(field_name{n}).day = PARAMS.ltsahd.day; % " "
    xwav_struct.(field_name{n}).hour = PARAMS.ltsahd.hour; % " "
    xwav_struct.(field_name{n}).minute = PARAMS.ltsahd.minute; % " "
    xwav_struct.(field_name{n}).second = PARAMS.ltsahd.secs; % " "
    xwav_struct.(field_name{n}).sample_freq = PARAMS.ltsahd.sample_rate; % xwav sampel frequency
    xwav_struct.(field_name{n}).file_name = PARAMS.ltsahd.fname; % xwav origin filename
    audio_data = audioinfo(files(n).name);
    xwav_struct.(field_name{n}).TotalSamples = audio_data.TotalSamples; % Samples per file
    xwav_struct.(field_name{n}).nxwav = PARAMS.ltsa.nrftot; % decimated xwavs per file
end


% Run Parameter setup
clear parm
[parm] = GPL_parameter_input_v3(PARAMS,xwav_struct);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear n






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Primary Loop - GPL Algorithm
% Ian Cosgrove

% Loop through all xwav files in se;ected folder.
for n = 1:length(files)
    
    % Loop through all disk writes within current xwav
    for q = 1:xwav_struct.(field_name{n}).nxwav
        
        % Display current window being processed
        fprintf('\n\nProcessing: %s\n',datestr(xwav_struct.(field_name{n}).julian_start_time(q)));

        % Determine sample range of current window
        SampleRangeStart = (parm.nrec*(q-1)) + 1;
        if q ~= xwav_struct.(field_name{n}).nxwav
            SampleRangeEnd = parm.nrec*(q) + 1;
        else
            SampleRangeEnd = parm.nrec*q; % Ensure requested samples ~> Total samples
        end
        
        % Extract audio data for the current window
        sub_data = audioread(files(n).name,[SampleRangeStart,SampleRangeEnd]);
        %%%%%%%%%%%%%%%%%%% audioread autonormalize: 'native' %%%%%%
        

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Error check for corrupted HARP data - Tyler Helble
        error_check = sort(abs(sub_data)); 
        error_check = mode(error_check);
        if error_check == 0
            fileID = fopen('errors.txt','at'); % Create text file
            fprintf(fileID, files(q).name,'\n'); % File Location
            fprintf(fileID, datestr(xwav_struct.(field_name{n}).julian_start_time(q)),'\n');
            fclose(fileID); % Close text file
            continue
        end
        
        % Kait: checking files have 0 data, updated function
        
        % End Helble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





        % Process current window through GPL algoirthm: identify calls,
        % extract call contours, take contour measurements, and perform
        % preliminary filtering.
        [GPL_struct] = GPL_v3(sub_data,parm);
        
        
        
    end % Individual xwav loop
    
end % Files loop
        


























% Set detection name
detection_name=(strcat('mod_01_08_24_detections_GPL_v2_',...
    num2str(parm.freq_lo),'_',num2str(parm.freq_hi),'_'));





toc
