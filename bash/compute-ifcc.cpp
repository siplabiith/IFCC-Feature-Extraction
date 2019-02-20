/*Brief Discription 

* This function extracts instantaneous frequency cepstral coefficients of speech data

* Usage: ifccExtract wavefile ifccfile

* Input Wavefile
* Output ifccfile

* Written By Dr.K.Srirama murthy
* Institute: Indian Institute of Technology, Hyderabad, India.
* Email: ksrm@iith.ac.in
* Date: 6 Nov 2012
*/


#include<iostream>
#include <fstream>
#include<sstream>
#include<stdlib.h>
#include<cmath>
#include<cstring>
#include<assert.h>
#include "fftw3.h"
#include <sndfile.h>

#define PI 3.141592653589793
#define BW 200.0
using namespace::std;

float *vector(int);
float **matrix(int,int);
void printVector(float *,int);
void printMatrix(float **,int,int);
void freeVector(float *);
void freeMatrix(float **,int);

float *wavread(char *, int, int *, int *);
void dither(float *,int,int);
void preemphasize(float *,int, float);

void melweights(int, int, int, float **,float *);
void linweights(int, int, int, float **,float *);
float mel2hz(float);
float hz2mel(float);
void smooth(float *,int, int);
float** spec2cep(float**, int, int,int);
float **deltas(float **, int,int,int);
void swap2(short *);
void swap4(int *);
void swapF(float *);
float *hamming(int);
void ifccExtract(char *, char *, int);
int main(int argc, char **argv)
{
	fftwf_init_threads();
	assert(("Usage: compute-ifcc wav.sp", argc==2));

        ifstream ifp;
        ifp.open(argv[1],ios::in);
	assert(("File does not exist",ifp));

        char file_id[1000];
        char buffer[1000];
        char file_name[1000];
        int chNum;
	string line;
	int count=0;	

        while(getline(ifp,line))
	{
               istringstream iss(line);
               if(!(iss>>file_id>>buffer>>buffer>>buffer>>buffer>>buffer>>chNum>>file_name>>buffer))
		{
			break;
		}
		ifccExtract(file_id, file_name, chNum-1);
		if(++count%10==0)
		{
			clog<<"LOG: Finished "<<count<<" files."<<endl;
		}
        }
}
void ifccExtract(char *file_id, char *file_name, int chNum)
{
        fftwf_plan_with_nthreads(1);

	float *wav;
	int numSamples=0;
	int fs=0;

	wav=wavread(file_name,chNum, &numSamples, &fs);
	if(wav==NULL)
	{
		return;
	}
	if(numSamples<500)
	{
		cerr<<"ERROR - Short duration file: "<<file_name<<endl;
		return;
	} 
	dither(wav,numSamples,1);
	preemphasize(wav,numSamples,0.97);
	

	int frameSize, frameShift, numFrames,L,numCepstra;
	frameSize=25*fs/1000;
	frameShift=10*fs/1000;
	numCepstra=20;
	numFrames=round(numSamples/frameShift);
	int numFramesFull=floor((numSamples-frameSize)/frameShift)+1;
	//numSamples=(numFrames-1)*frameShift+frameSize;

// COMPUTE DFT OF THE GIVEN SPEECH SIGNAL.	
	fftwf_complex *X;
	fftwf_plan p;
	L=floor(numSamples/2)+1;
	X = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*L);
	p= fftwf_plan_dft_r2c_1d(numSamples, wav, X, FFTW_ESTIMATE );
	fftwf_execute(p);
	fftwf_destroy_plan(p);
	freeVector(wav);
	
// CONSIDER THE BAND OF INTREST. 	

	float minFreq=20.0;
	float maxFreq=3850.0;
        int minIdx, maxIdx;
	minIdx=floor(2*L*minFreq/fs);
	maxIdx=floor(2*L*maxFreq/fs);
	L=(maxIdx-minIdx);
//      cout<<numSamples<<" "<<L<<endl<<endl;

//COMPUTE BAND WEIGHTS

	int numBands=40; //ceil(hz2mel(fs/2))+1;
        float **W=matrix(numBands,L);
        float *midFreq=vector(numBands);
//	melweights(L,fs,numBands,W, midFreq);
        linweights(L,fs,numBands,W,midFreq);
        
        //printVector(midFreq,numBands);
	fftwf_complex *Y, *ansig, *dsig;
	Y = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*numSamples);
	ansig=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*numSamples);
	dsig=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*numSamples);

	fftwf_plan p1,p2;
	p1=fftwf_plan_dft_1d(numSamples, Y, ansig,FFTW_BACKWARD,FFTW_ESTIMATE );	
	p2=fftwf_plan_dft_1d(numSamples, Y, dsig,FFTW_BACKWARD,FFTW_ESTIMATE );	

	float *hamwin=hamming(frameSize);
	float **ifSpec=matrix(numBands, numFrames);
	float *instFreq=vector(numSamples);
	float *hilenv=vector(numSamples);
	float **cepstra;
        for(int i=0;i<numSamples;i++)
	{
	    Y[i][0]=0;
	    Y[i][1]=0;
	}
	for(int i=0;i<numBands;i++)
	{
		midFreq[i]=midFreq[i]+minIdx;
		for(int j=minIdx;j<maxIdx;j++)
		{
                        Y[j][0]=X[j][0]*W[i][j-minIdx];
			Y[j][1]=X[j][1]*W[i][j-minIdx];
		}
		fftwf_execute(p1);
		
                for(int j=minIdx;j<maxIdx;j++)
		{
			Y[j][0]*=j;
			Y[j][1]*=j;
		}
		fftwf_execute(p2);
		
		for(int j=0; j<numSamples;j++)
		{
			instFreq[j]=ansig[j][0]*dsig[j][0]+ansig[j][1]*dsig[j][1];
			hilenv[j]=ansig[j][0]*ansig[j][0]+ansig[j][1]*ansig[j][1];
		}
		smooth(instFreq, numSamples, frameSize+1);
		smooth(hilenv, numSamples, frameSize+1);
		for(int j=0;j<numSamples;j++)
		{
			instFreq[j]/=hilenv[j];
                        instFreq[j]-=midFreq[i];
			instFreq[j]*=(2.0*PI)/numSamples;
		}
                //printVector(instFreq,numSamples);
		
		
		for(int j=0;j<numFramesFull; j++)
		{
			ifSpec[i][j]=0;
			float wtNorm=0;
			for(int k=0; k<frameSize; k++)
			{
				ifSpec[i][j]+=instFreq[j*frameShift+k]*hamwin[k];
			}
			ifSpec[i][j]*=100.0;
		}
		for(int j=numFramesFull;j<numFrames; j++)
		{
			ifSpec[i][j]=0;
			float wtNorm=0;
			for(int k=0; k<frameSize&&j*frameShift+k<numSamples; k++)
			{
				ifSpec[i][j]+=instFreq[j*frameShift+k]*hamwin[k];
				wtNorm+=hamwin[k];	
			}
			ifSpec[i][j]=ifSpec[i][j]*100.0/wtNorm;
		}

	}
	fftwf_destroy_plan(p1);
	fftwf_destroy_plan(p2);
	fftwf_free(X);
	fftwf_free(Y);
	fftwf_free(ansig);
	fftwf_free(dsig);
	freeVector(instFreq);
	freeVector(hilenv);
	freeMatrix(W,numBands);
	freeVector(hamwin);
	
	cepstra=spec2cep(ifSpec,numBands, numFrames,numCepstra);
	cout<<file_id<<" ["<<endl;
        printMatrix(cepstra,numFrames,numCepstra);
	cout<<"]"<<endl;
	freeMatrix(ifSpec,numBands);
	freeMatrix(cepstra,numFrames);
}
float *hamming(int N)
{
    float alpha=0.54;	
    int half; 
    float temp;
    float *x; 
    x=vector(N);	 
    for (int i =0; i<N/2.0; i++)
    {
	x[i] = alpha - (1 - alpha)*cos(2*M_PI*i/(N-1));
	x[N-1-i]=x[i];
    }	
    float norm=0;	
    
    for(int i=0; i<N; i++)
    {
        norm+=x[i];
    }	

    for(int i=0; i<N; i++)
    {
        x[i]/=norm;
    } 
	
    return x;
}
void swap2(short *p) 
{
	char temp,*q;
   	q = (char*) p;
   	temp = *q; *q = *(q+1); *(q+1) = temp;
}

void swap4(int *p)
{
	char temp,*q;
	q = (char*) p;
   	temp = *q; *q = *(q+3); *(q+3) = temp;
   	temp = *(q+1); *(q+1) = *(q+2); *(q+2) = temp;
}
void swapF(float *p)
{
	char temp,*q;
	q = (char*) p;
   	temp = *q; *q = *(q+3); *(q+3) = temp;
   	temp = *(q+1); *(q+1) = *(q+2); *(q+2) = temp;
}

float **deltas(float **x, int numFrames, int numCepstra, int winLen)
{

	float **d=matrix(numFrames,numCepstra);
    	
        int hlen = floor(winLen/2);
        float normFactor=0;
        for(int i=-hlen;i<=hlen;i++)
        {
            normFactor+=i*i;
        }
	
        for(int i=0;i<numFrames;i++)
	{
		for(int j=0; j<numCepstra;j++)
		{
			d[i][j]=0;
			for(int k=-hlen;k<=hlen;k++)
			{
                            if(i+k<0)         
                                d[i][j]+=x[0][j]*k;
                            else if (i+k>=numFrames)
                                d[i][j]+=x[numFrames-1][j]*k;
                            else
                                d[i][j]+=x[i+k][j]*k;
			}
                        d[i][j]/=normFactor;
		}
	}			
	return (d);
}


float **spec2cep(float **spec, int numBands, int numFrames, int numCepstra)
{
	float **dctm;
	float **cepstra;
	dctm=matrix(numCepstra, numBands);
	cepstra=matrix(numFrames,numCepstra);

	for ( int i = 0; i < numCepstra; i++ )
        {
        	for ( int j = 0; j < numBands; j++ )
		{
                	dctm[i][j] = cos(PI*i*(2.0*(float)j+1)/(2.0*(float)numBands)) * sqrt(2.0/(float)numBands);
        	}
	}
	
	for(int j=0; j<numBands; j++)
	{
		dctm[0][j]/=sqrt(2);
	}

	for(int i=0; i<numFrames; i++)
	{
		for(int j=0; j<numCepstra; j++)
		{
			cepstra[i][j]=0;
			for(int k=0;k<numBands;k++)
			{
				//cepstra[i][j]+=log(spec[k][i])*dctm[j][k];
				cepstra[i][j]+=(spec[k][i])*dctm[j][k];
			}
			//cepstra[i][j]*=pow(j+1,0.6);	
		}
	}	
	return(cepstra);	
}
void smooth(float *sig, int numSamples, int winSize)
{
	int d;
	d=floor(winSize/2);
	winSize=2*d+1;
	float *x;
	x=vector(numSamples+winSize-1);
	for(int i=0; i<d;i++)
	{
		x[i]=0;
		x[numSamples+i]=0;
	}
	for(int i=0;i<numSamples;i++)
	{
		x[i+d]=sig[i];
		sig[i]=0.0;
	}
	for(int i=d; i<numSamples+d;i++)
	{
		for(int j=i-d; j<=i+d; j++)
		{
			sig[i-d]+=x[j];
		}
	}

	freeVector(x);
}	

void melweights(int L, int fs, int numBands, float **W, float *midFreq)
{
	float nyqmel;
       	nyqmel= hz2mel(fs/2);
    	float step_mels;
        step_mels = nyqmel/(numBands - 1);
    	float *binmels;
       	binmels=vector(L);
	float dB=48.0;

    	for ( int i = 0; i < L; i++ )
    	{
		binmels[i] = hz2mel(((float)i*(fs/2.0))/(L-1));
    	}

    	for ( int i = 0; i < numBands; i++ )
    	{
		float f_mel_mid = i*step_mels;
		for ( int j = 0; j < L; j++ )
		{
	    		W[i][j] = exp(-0.5*pow(binmels[j]-f_mel_mid,2));
		}
                midFreq[i]=round(mel2hz(f_mel_mid)/fs*L*2.0);
    	}

    	freeVector(binmels);
}
void linweights(int L, int fs, int numBands, float **W,float *midFreq)
{
    	float step = 1.0/(numBands - 1);
    	float *bin;
	float sigma_sqr=pow(BW/(fs*0.8325546),2);
       	bin=vector(L);
    	for ( int i = 0; i < L; i++ )
    	{
		bin[i] = (float) i/(L-1);
    	}

    	for ( int i = 0; i < numBands; i++ )
    	{
		float f_mid = i*step;
		for ( int j = 0; j < L; j++ )
		{
	    		W[i][j] = exp(-pow(bin[j]-f_mid,2)/(2*sigma_sqr));
		}
                midFreq[i]=round(i*step*L);
    	}

    	freeVector(bin);
}


float hz2mel( float hz )
{
	float f_0 = 0; // 133.33333;
   	float f_sp = 200.0/3.0; // 66.66667;
    	float brkfrq = 1000;
    	float brkpt  = (brkfrq - f_0)/f_sp;  
    	float z;  
    	float logstep = exp(log(6.4)/27.0);

    	// fill in parts separately
    	if (hz < brkfrq)
    	{
		z = (hz - f_0)/f_sp;
    	}
    	else
    	{
		z  = brkpt+((log(hz/brkfrq))/log(logstep));
    	}
    	return z; 
}

float mel2hz(float mel)
{
	float f_0 = 0;
  	float f_sp = 200.0 / 3.0;
  	float brkfrq = 1000;
  	float brkpt = (brkfrq - f_0)/f_sp;
  	float z;
  	float logstep = exp(log(6.4)/27);

  	if (mel < brkpt)
  	{
      		z = f_0 + f_sp * mel;
  	}
  	else
  	{
      		z = brkfrq * exp(log(logstep) * (mel - brkpt));
  	}
  	return z;
}


void preemphasize(float *wav, int numSamples, float alpha)
{
	if(alpha<0|alpha>1.0)
	{
		cerr<<"Preemphasis coeffient must be in the range 0 to 1"<<endl;
		return;
	}

	for(int i=numSamples-1; i>0;i--)
	{
		wav[i]=wav[i]-alpha*wav[i-1];
	}
}




void dither(float *wav,int numSamples,int scale)
{
	float r;
	for(int i=0;i<numSamples;i++)
	{
		r=((float) rand())/RAND_MAX;
		wav[i] +=round(scale*(2*r-1));
	}
}
float *wavread(char *wavFile, int chNum, int *numSamples, int *fs)
{
        SNDFILE *sf;
        SF_INFO info;
        sf = sf_open(wavFile,SFM_READ,&info);
        if (sf == NULL)
        {
                cerr<<"ERROR - Failed to open the file:"<<wavFile<<endl;
                return(NULL);
        }
        int num_items=info.channels*info.frames;
        short *buffer=new short [num_items];
        sf_read_short(sf,buffer,num_items);
        float *wav=vector(info.frames);
        for(int i=chNum, j=0; i<num_items; i+=info.channels, j++)
        {
                wav[j]=buffer[i];
        }
        delete [] buffer;
        *fs=info.samplerate;
	*numSamples=info.frames;
	return(wav);
}
float *vector(int length)
{
		float *vec;
		vec=(float *)malloc(length*sizeof(float));
		if(!vec)
		{
				printf("Memory allocation failed in vector().\n");
				exit(1);
		}
		
		return(vec);
}

float **matrix(int row, int col)
{
		float **mat;
		//mat=new float* [row];
		mat=(float **)malloc(row*sizeof(float *));
		if(!mat)
		{
				printf("Memory allocation failed in Matrix().\n");
				exit(1);
		}
		
		for(int i=0;i<row;i++)
		{
				//mat[i]=new float [col];
				mat[i]=(float *)malloc(col*sizeof(float));
				if(!mat[i])
				{
					printf("Memory allocation failed in Matrix().\n");
					exit(1);
				}
		}
		return(mat);
}

void freeVector(float *vec)
{
		free(vec);
}

void freeMatrix(float **mat, int row)
{
		for(int i=0;i<row;i++)
		{
				free(mat[i]);
		}
		free(mat);
}
void printVector(float *vec, int length)
{
		for(int i=0;i<length;i++)
		{
				cout<<vec[i]<<" ";
		}
		cout<<endl;
}

void printMatrix(float **mat, int row, int col)
{
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			cout<<mat[i][j]<<" ";
		}
			cout<<endl;
	}
}
