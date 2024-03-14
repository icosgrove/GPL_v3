GPL_v3 is a reconstructed and documented version of Tyler Helble's Generalized Power-Law (GPL) Detector Algorithm.
The theory behind the detector can be read in Helble et al. "A generalized power-law detection 
algorithm for humpback whale vocalizations" (2012) at https://www.cetus.ucsd.edu/publications.html

This is a work-in-progress and may not run/produce desired outputs in it's current state. (As of 03/13/2024)

Algorithm Summary: 

The GPL detector is designed to detect tonal underwater marine mammal acoustic signals from HARP data in environments that vary in background noise level. The detector will parse through x.wav files in short time windows and use the power-law algorithm 
to locate regions of time with energy spikes over a specified frequency range. The detector has functions to remove noise from spectrograms and process only portions of them that are considered 'signal present'. The algorithm uses
a noise floor value calculated from the mean background noise level of a given window of data to determine if an energy spike is strong enough to be a marine mammal vocalization. This noise floor is adjusted every time a new window 
is processed to accurately reflect the changing underwater environment. 

The detector will produce regions of time within a given audio file that have enough energy relative to the background noise to be a potential marine mammal vocalization. Then, the regions of time will be analyzed further to determine
where in frequency the energy came from, which produces a contour for a given detection. These contours are measured and filtered based on specifications set by the user to look for a specific signal type. 

Outputs include the detection start/end time, associated contour, and relevant contour/background noise measurements.
