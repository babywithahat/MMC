//MASTER CHANNEL
#include <stdio.h>
#include <stdlib.h>
#include <sf.h>
#include <sfmisc.h>

#include "fm.h"
#include "cat.h"
//#include "file2array.h"
#include "mix.h"
#include "stripchannel.h"
//#include "duroftone.h"
//#include "allocsetup.h"
//#include "loadarray.h"

//call module functions
int main(void){
	
	int i;
	//declare channel count and sampling rate
	int nchan = 8,
		sr = 44100;
	
	//declare master db
	double dB = 90;
	
	//REVISE. declare tone duration, start time, and total dur
	double totalDur=0; //REVISE

	//declare master output.
	short *output;
	
	//start time of tones
	double startTime[] = {0, 0, 15.7, 27.7, 27.7, 27.7, 31.7, 35.0, 38.3, 38.3, 90, 51.5, 7.7, 123.36075, 23.7, 123.36075, 123.36075, 123.76075};
	int nTones = 18;
	
	//allocate memory for output.
	double toneDur;
	for(i = 0; i < nTones; i++){
		printf("Running duroftone()...\n");
		toneDur = durOfTone(i+1, startTime[i]);
		if(toneDur >= totalDur)
			totalDur = toneDur;
		
		printf("Done. \n/////// totalDur = %f\n", totalDur);
	}
	totalDur = totalDur + 1;	
	printf("Allocating Memory for output...");
	output = (short *)Malloc(totalDur * sr * nchan * sizeof(short));	
	
	printf("Done.\n");
	
	//allocate mix output
	
	printf("Running mixSetup()...");
	
	mixSetup(totalDur, sr, nchan);
	
	printf("Done.\n");
	
	//lets do this thing
	
	printf("Declaring malloc multipliers...");
	
	int mallocCar, mallocMod, mallocPD, mallocD, mallocAmp, mallocT,
		mallocA, mallocPan;
	
	printf("Done.\n");
	
	double *Carrier;
	double *Modulator;
	double *PD;
	double *D;
	double *Amp;
	double *T;
	double *A;
	double *xpan;
	double *ypan;
	for(i=0; i<nTones; i++){
		printf("///////\nTONE %d\n///////\n",i+1);
		printf("Running allocSetup()...\n");	
	
		allocSetup(i+1, &mallocCar, &mallocMod, &mallocPD, &mallocD,
			&mallocAmp, &mallocT, &mallocA, &mallocPan);
	
		printf("Done.\n");
		printf("mallocPan = %d\n", mallocPan);
		printf("Allocating Memory for double arrays...");
	
		Carrier = (double *)Malloc(mallocCar * sizeof(double));
		Modulator = (double *)Malloc(mallocMod * sizeof(double));
		PD = (double *)Malloc(mallocPD * sizeof(double));
		D = (double *)Malloc(mallocD * sizeof(double));
		Amp = (double *)Malloc(mallocAmp * sizeof(double));
		T = (double *)Malloc(mallocT * sizeof(double)); 
		A = (double *)Malloc(mallocA * sizeof(double));
	
		printf("mallocPan = %d", mallocPan);
	
		xpan = (double *)Malloc(mallocPan * sizeof(double));
		ypan = (double *)Malloc(mallocPan * sizeof(double));
	
		printf("Done.\n");
	
		int nFMTones;
		int	nPoints;
	
		printf("Running loadarray()...\n");
	
		loadarray(i+1, Carrier,Modulator,PD,D,Amp,&nFMTones,
			T, A, &nPoints, xpan, ypan);
	
		printf("Done.\n");
		printf("Running stripchannel()...\n");
	
		stripchannel(dB,startTime[i],Carrier,Modulator,PD,D,Amp,&nFMTones,
			T,A,&nPoints,0,xpan, ypan);
	
		printf("Done.\n");
	
		Free(Carrier);
		Free(Modulator);
		Free(PD);
		Free(D);
		Free(Amp);
		Free(T);
		Free(A);
		Free(xpan);
		Free(ypan);
	}
	// convert the mix output from doubles to shorts
	output = mixConvertTo16Bits();
	
	//save the tones
	sfsave("week5.wav", output, totalDur, sr, nchan);

	//release memory
	Free(output);

	exit(EXIT_SUCCESS);
}
