For extracting IFCC features run below commands


data/train ---> data directory of wavfiles

make_ifcc  ----> log directory for checking the progress of feature extraction

ifcc --------> output directory where the IFCC features will store

you can change number of jobs according to your system


python IFCC_features_extract.py <data> <log> <output> <num_jobs>

example:

python IFCC_features_extract.py data/train make_ifcc ifcc 1

NOTE:-

if you want to change output format, you can change output file format in IFCC_features_run.py

wav.scp:
	uttID channel wave_file_location
