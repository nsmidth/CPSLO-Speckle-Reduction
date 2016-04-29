all:
	gcc -o fftw_psd.so -shared -fPIC fftw_psd.c
	gcc -c -I/$HOME/include fftw_psd.c
	gcc fftw_psd.o -lfftw3 -lm
	rm fftw_psd.o
	mv a.out fftw_psd
