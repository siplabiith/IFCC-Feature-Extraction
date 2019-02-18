# IFCC-Feature-Extraction

If you use this code 

please cite us

#References

<a href="https://www.researchgate.net/publication/283203401_Epoch_Extraction_by_Phase_Modelling_of_Speech_Signals">1. Karthika Vijayan and K Sri Rama Murty, "Epoch Extraction by Phase Modelling of Speech Signals," Circuits, Systems, and Signal Processing, vol. 35, no. 7, pp. 2584-2609, July 2016</a>
  
<a href="https://ieeexplore.ieee.org/document/8308665">2. Shekhar Nayak, Saurabhchand Bhati and K Sri Rama Murty, "An Investigation into Instantaneous Frequency Estimation Methods for Improved Speech Recognition Features," in Proc. IEEE Global Conference on Signal and Information Processing (GlobalSIP), Nov. 14-16, 2017, Montreal, Canada.</a>

<a href="https://www.semanticscholar.org/paper/Feature-extraction-from-analytic-phase-of-speech-Vijayan-Kumar/41f56c120d0657c1ba500b1d39ad0d856b152dbb">3. Karthika Vijayan, Vinay Kumar and K. Sri Rama Murty, “Feature extraction from analytic phase of speech signals for speaker verification”, in INTERSPEECH, September 2014, Singapore, pp.1658-1662.</a>

We are providing IFCC feature extraction in 2 languages

1.python<br></br
2.bash

Follow these instructions to extract IFCC features

##using kaldi

Libraries to install

	1.fftw3
	
	2.libsnd
	
after installing the library files comile the code ( ./compile.sh)

for extracting IFCC features for .sph files

	comput-ifcc.cpp
	
for extracting IFCC features for wav files

	compute-ifcc_wav.scp

use make_ifcc.sh code to extract IFCC features

change the comput-ifcc.cpp file for changing the number of coefficients

##using python

python IFCC_features_extract.py <data> <log> <output> <num_jobs>
  
example:

where

data       ----> data directory of wavfiles

make_ifcc  ----> log directory for checking the progress of feature extraction

ifcc       ----> output directory where the IFCC features will store

num_jobs   ----> no of jobs

python IFCC_features_extract.py data/train make_ifcc ifcc 1

NOTE:-

if you want to change output format, you can change output file format in IFCC_features_run.py

you can change number of jobs according to your system

wav.scp:

	uttID channel wave_file_location
