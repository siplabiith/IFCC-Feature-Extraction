#code for extracting IFCC features
# importing the modules 
import datetime
import sys
from scipy.io import wavfile
import numpy as np
import logging
import random
import math
import threading
import multiprocessing
from IFCC_features_fun import *
#using multiprocessing to obtain IFCC features
base_name=sys.argv[1]  #taking data dir from terminal
log_dir=sys.argv[2]    #taking logdir from terminal
ifcc_dir=sys.argv[3]   #taking ifccdir from terminal
number_jobs=int(sys.argv[4]) #taking number of jobs from terminal
#getting the wav.sp files for different jobs
file_names=[]



for i in range(1,int(number_jobs)+1,1):
	name=log_dir+"/"+"wav_"+base_name+".%s.scp"%i
	file_names.append(name)


location=base_name+"/wav.scp"
lines=open(location).read().splitlines()
lines=np.array(lines)
total_lines=float(len(lines))
files_job=int(math.ceil(total_lines/number_jobs))
for i in range(0,number_jobs,1):
	lines1=lines[i*files_job:(i+1)*files_job]
	name=log_dir+"/"+"wav_"+base_name+".%s.scp"%(i+1)
	np.savetxt(name,lines1,fmt="%s")


#getting the log files for different jobs
log_names=[]
for i in range(1,int(number_jobs)+1,1):
	name=log_dir+"/make_ifcc.%s.log"%i
	log_names.append(name)


#getting output file names for storing ifcc features
ifcc_names=[]
for i in range(1,int(number_jobs)+1,1):
	name=ifcc_dir+"/ifcc.%s.features"%i
	ifcc_names.append(name)


process=[None]*number_jobs
for i in range(0,number_jobs,1):
	process[i]=multiprocessing.Process(target=IFCC_features,args=(file_names[i],ifcc_names[i],log_names[i],))
for i in range (0,number_jobs,1):
	process[i].start()
for i in range (0,number_jobs,1):
	process[i].join()
