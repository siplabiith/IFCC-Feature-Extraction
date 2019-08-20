# Copyright Speech Information Processing Lab, IIT Hyderabad .  All Rights Reserved.

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
PI=3.141592653589793
M_PI=3.141592653589793
BW=200.0
b=1
def movingaverage(values,window):
	weights=np.repeat(1.0,window)/window
	sma=np.convolve(values,weights,'same')
	return sma
def mel2hz(mel):
	f_0=0
	f_sp=200.0/3.0
	brkfrq=1000
	brkpt=(brkfrq-f_0)/f_sp
	logstep=math.exp(math.log(6.4)/27)
	if(mel<brkpt):
		z=f_0+f_sp*mel
	else:
		z=brkfrq*math.exp(math.log(logstep)*(mel-brkpt))
	return z
def hz2mel(hz):
	f_0=0
	f_sp=200.0/3.0
	brkfrq=1000
	brkpt=(brkfq-f_0)/f_sp
	logstep=np.exp(np.log(6.4)/27.0)
	if(hz<brkfrq):
		z=(hz-f_0)/f_sp
	else:
		z=brkpt+((np.log(hz/brkfrq))/np.log(logstep))
	return z
def dither(wav,numSamples,scale):
	for i in range(0,numSamples,1):
		r=random.random()
		wav[i]=wav[i]+round(scale*(2*r-1))
	return wav


def preemphasis(wav,numSamples,alpha):
	if alpha < 0 or alpha > 1.0 :
		print ( "ERROR :- Preemphasis coeffient must be in the range 0 to 1 ")
		exit()
	for i in range(numSamples-1,0,-1):
		wav[i]=wav[i]-alpha*wav[i-1]
	return wav
def hamming(N):
	alpha=0.54
	x=vector(N)
	for i in range(0,int(N/2.0),1): 
		x[i]=alpha-(1-alpha)*math.cos(2*M_PI*i/(N-1))
		x[N-1-i]=x[i]
	norm=0.0
	for i in range(0,N,1):
		norm=norm+x[i]
	for i in range(0,N,1):
		x[i]=x[i]/norm
	return x
def matrix(a,b):
	matrix1=[]
	for i in range(0,a,1):
		matrix1.append(vector(b))
	return matrix1
def vector(a):
	vec=[]	
	for i in range(0,a,1):
		vec.append(0.0)
	return vec

def linweights(L,fs,numBands,W,midFreq):
	W=np.array(W)
	step=float(1.0/(numBands-1))
	bin1=np.zeros(L)
	sigma_sqr=float(pow(BW/(fs*0.8325546),2))
	for i in range(0,L,1):
		bin1[i]=float(i/(L-1.0))
	for i in range(0,numBands,1):
		f_mid=float(i*step)
		W[i,:]=np.exp(((bin1-f_mid)**2)/(-2.0*sigma_sqr))
		midFreq[i]=round(i*step*L)
	return W,midFreq

def smooth(sig,numSamples,winSize):
	d=math.floor(winSize/2)
	winSize=2*d+1
	sig=movingaverage(sig,winSize)
	return sig

def spec2cep(spec,numBands,numFrames,numCepstra):
	dctm=matrix(numCepstra,numBands)
	cepstra=matrix(int(numFrames),numCepstra)
	for i in range(0,numCepstra,1):
		for j in range(0,numBands,1):
			dctm[i][j]=math.cos(PI*i*(2.0*float(j)+1)/(2.0*float(numBands)))*math.sqrt(2.0/float(numBands))	
	for j in range(0,numBands,1):
		dctm[0][j]=dctm[0][j]/math.sqrt(2)
	for i in range(0,int(numFrames),1):
		for j in range(0,numCepstra,1):
			cepstra[i][j]=0
			for k in range(0,numBands,1):
				cepstra[i][j]=cepstra[i][j]+(spec[k][i])*dctm[j][k]
	return cepstra
def melweights(L,fs,numBands,W,float,midFreq):
	nyqmel=hz2mel(fs/2)
	step_mels=nyqmel/(numBands-1)
	binmels=vector(L)
	dB=48.0
	for i in range(0,L,1):
		binmels[i]=hz2mel((float(i)*(fs/2.0))/(L-1))
	for i in range(0,numBands,1):
		f_mel_mid=i*step_mels
		for j in range(0,L,1):
			W[i][j]=np.exp(-0,5*pow(binmels[j]-f_mel_mid,2))
		midFreq[i]=round(mel2hz(f_mel_mid)/fs*L*2.0)
	return W,midFreq
def deltas(x,numFrames,numCepstra,winLen):
	d=matrix(numFrames,numCepstra)
	hlen=math.floor(winLen/2)
	normFactor=0
	for i in range(-hlen,hlen+1,1):
		normFactor=normFactor+i*i
	for i in range(0,numFrames,1):
		for j in range(0,numCepstra,1):
			d[i][j]=0
			for k in range(-hlen,hlen+1,1):
				if(i+k<0):
					d[i][j]=d[i][j]+x[0][j]*k
				elif (i+k>=numFrames):
					d[i][j]=d[i][j]+x[numFrames-1][j]*k		
				else:
					d[i][j]=d[i][j]+x[i+k]*k
			d[i][j]=d[i][j]/normFactor
	return d
def narrowBand(X,W,midFreq,hamwin,Y,Y1):
	global ifSpec
	ifSpec=matrix(numBands,int(numFrames))
	for i in range(0,numBands,1):
			midFreq[i]=midFreq[i]+minIdx
			for j in range(int(minIdx),int(maxIdx),1):
				Y[j]=X[j]*W[i][j-int(minIdx)]
			for j in range(int(minIdx),int(maxIdx),1):
				Y1[j]=j*Y[j]
			Y=np.array(Y)
			Y1=np.array(Y1)
			ansig1=np.fft.ifft(Y,numSamples,axis=-1)
			dsig1=np.fft.ifft(Y1,numSamples,axis=-1)
			ansig=[ansig1.real,ansig1.imag]
			dsig=[dsig1.real,dsig1.imag]
			instFreq=(numSamples*ansig[0])*(numSamples*dsig[0])+(numSamples*ansig[1])*(numSamples*dsig[1])
			hilenv=(numSamples*ansig[0])*(numSamples*ansig[0])+(numSamples*ansig[1])*(numSamples*ansig[1])
			instFreq=smooth(instFreq,numSamples,frameSize+1)
			hilenv=smooth(hilenv,numSamples,frameSize+1)
			instFreq=instFreq//hilenv
			instFreq=instFreq-midFreq[i]
			instFreq=instFreq*((2.0*PI)/numSamples)
			for j in range(0,int(numFramesFull),1):
				instFreq1=instFreq[j*frameShift:(j*frameShift)+(frameSize):1]
				ifSpec[i][j]=np.matmul(instFreq[j*frameShift:(j*frameShift)+(frameSize):1],hamwin)
				ifSpec[i][j]=ifSpec[i][j]*100.0
			for j in range(int(numFramesFull),int(numFrames),1):
				ifSpec[i][j]=0
				wtNorm=0.0
				k=0
				while k<frameSize and j*frameShift+k<numSamples :
					ifSpec[i][j]=ifSpec[i][j]+instFreq[j*frameShift+k]*hamwin[k]
					wtNorm=wtNorm+hamwin[k]
					k=k+1
				ifSpec[i][j]=(ifSpec[i][j]*100.0)/wtNorm
	del(instFreq)
	del(hilenv)
	del(ansig1)
	del(dsig1)
	del(Y1)
	del(Y)
	del(ansig)
	del(dsig)
	return ifSpec
#print(len(sys.argv))

#it will take three arguments which are given by main script
numBands=40
minFreq=20.0
maxFreq=3850.0
numCepstra=20
def IFCC_features(file_name,output_file,log_file):
	erase=open(log_file,"w")
	logging.basicConfig(filename=log_file,level=logging.DEBUG)
	logger = logging.getLogger()
	f=open(output_file,'w')
	f=open(output_file,'ab')
	NEWLINE_SIZE_IN_BYTES = -1
	file = open(file_name)
	for line in file:
		words=line.split(" ")
		file_id=words[0]
		numChan=int(words[1])
		inter=words[2].split("\n")
		file_name=inter[0]
		if ".sph" in file_name:
			new_file_name=utter_path.replace("sph","wav")
			cmd="sox %s %s"%(file_name,new_file_name)
			a=os.system(cmd)
			fs, data1 = wavfile.read(new_file_name)
			if isinstance(data1[0],np.ndarray):
				data=data1[:,numChan]
			else:
				data=data1
			os.remove("%s"%new_file_name)
		else:
			fs, data1 = wavfile.read(file_name)
			if isinstance(data1[0],np.ndarray):
				data = data1[:,numChan]
			else:
				data=data1
		wav=data.tolist()
		global numSamples
		numSamples=len(wav)
		if(len(data) < 500):
			print ("WARNING CHECK examples.log file")
			logger.warning('number of samples in the audion file are less than 500')
			exit()
		#wav=dither(wav,numSamples,1)
		wav=preemphasis(wav,numSamples,0.97)
		global frameSize
		frameSize=25*fs/1000 #25 ms
		global frameShift
		frameShift=10*fs/1000 #10 ms
		global numFrames
		global numFramesFull
		numFrames=round(numSamples/frameShift)  # this number may not include last frame 
		numFramesFull=math.floor((numSamples-frameSize)/frameShift)+1 #this number will definetly include last frame
		L=math.floor(numSamples/2)+1
		X=np.fft.fft(wav,numSamples)
		global minIdx
		minIdx=math.floor(2*L*minFreq/fs)
		global maxIdx
		maxIdx=math.floor(2*L*maxFreq/fs)
		L=(maxIdx-minIdx)
		W=matrix(numBands,int(L))
		midFreq=vector(numBands)
		W,midFreq=linweights(int(L),fs,numBands,W,midFreq)
		instFreq=vector(numSamples)
		hilenv=vector(numSamples)
		hamwin=np.hamming(frameSize)
		hamwin=hamwin/sum(hamwin)
		Y=vector(numSamples)
		Y1=vector(numSamples)
		ifSpec=narrowBand(X,W,midFreq,hamwin,Y,Y1)
		cepstra=spec2cep(ifSpec,numBands, numFrames,numCepstra)
		f.write("%s [ "%file_id)
		np.savetxt(f,cepstra,fmt='%.4f')
		f.seek(NEWLINE_SIZE_IN_BYTES, 2)
		f.truncate()
		f.write(" ] ")
		f.write("\n")
		print(cepstra)
		logger.info("processed IFCC features successfully for %s "%file_id)
