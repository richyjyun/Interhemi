# [Interhemi](https://www.frontiersin.org/articles/10.3389/fnins.2021.782188/full)

Code used for analyzing data for ["Cortical stimulation paired with volitional unimanual movement affects interhemispheric communication."](https://www.frontiersin.org/articles/10.3389/fnins.2021.782188/full) 

Data consisted of neural data from 32 channels of a epidural implant collected by g.USBAMP (g.tec), and wrist acceleration data collected through a NI-DAQ device using custom MATLAB code.

+u contains mainly utility functions including loading the metadata from a Google spreadsheet shared between investigators, calculating basic metrics such as reaction time, aligning behavior data to the neural data, and compiling a session list with all experiments.

+a contains scripts used for analzying and plotting the LFP analyses. Granger's causality was obtained using the MVGC toolbox (Barnett & Seth, 2014) and power over time was calculated in part with the Chronux package (http://chronux.org/; Mitra & Bokil, 2008) in conjunction with custom code.


[1] L. Barnett and A. K. Seth, The MVGC Multivariate Granger Causality Toolbox: A New Approach to Granger-causal Inference, J. Neurosci. Methods 223, 2014.

[2] P. Mitra and H. Bokil, Observed Brain Dynamics. Oxford University Press, 2008

Figure demonstrating intra- and interhemispheric interactions during unimanaul movement. 

<p align="center">
  <img width="750" height="302.5" src="https://github.com/richyyun/Interhemi/blob/main/Body%20-%20MovementModel.png">
</p>
