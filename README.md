# Interhemi

Part of code used for analyzing data for "Cortical stimulation paired with volitional unimanual movement affects interhemispheric communication." Manuscript has been accepted to Frontiers in Neuroscience.

Data consisted of neural data from 32 channels of a epidural implant collected by g.USBAMP (g.tec), and wrist acceleration data collected through a NI-DAQ device using custom MATLAB code.

+u contains mainly utility functions including loading the metadata from a Google spreadsheet shared between investigators, calculating basic metrics such as reaction time, aligning behavior data to the neural data, and compiling a session list with all experiments.

+a contains scripts used for analzying and plotting the LFP analyses. Granger's causality was obtained using the MVGC toolbox (Barnett & Seth, 2014) and power over time was calculated in part with the Chronux package (http://chronux.org/; Mitra & Bokil, 2008) in conjunction with custom code.


L. Barnett and A. K. Seth, The MVGC Multivariate Granger Causality Toolbox: A New Approach to Granger-causal Inference, J. Neurosci. Methods 223, 2014.

P. Mitra and H. Bokil, Observed Brain Dynamics. Oxford University Press, 2008

<p align="center">
  <img width="500" height="300" src="https://github.com/richyyun/Interhemi/blob/main/Figure%201%20-%20MovementModel.tif">
</p>
