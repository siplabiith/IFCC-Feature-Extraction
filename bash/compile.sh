g++ -openmp compute-ifcc-test.cpp  -O3 -o compute-ifcc-test  -I/home/rafi/kaldi/fftw3/include -I/home/rafi/kaldi/libsnd/include -L/home/rafi/kaldi/fftw3/lib/ -L/home/rafi/kaldi/libsnd/lib -lfftw3f -lfftw3f_threads -lm -lsndfile -lpthread