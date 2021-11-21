# Interhemi

Part of code used for analyzing data for "Cortical stimulation paired with volitional unimanual movement affects interhemispheric communication." Manuscript has been accepted to Frontiers in Neuroscience.

Data consisted of neural data from 32 channels of a epidural implant collected by g.USBAMP (g.tec), and wrist acceleration data collected through a NI-DAQ device using custom MATLAB code.

+u contains mainly utility functions including loading the metadata from a Google spreadsheet shared between investigators, calculating basic metrics such as reaction time, aligning behavior data to the neural data, and compiling a session list with all experiments.

+a contains scripts used for analzying and plotting the LFP analyses. 
